# Estimators - DTFE density and velocity estimation functions

"""
    AbstractEstimator{T}

Abstract base type for DTFE estimators.
Subtypes implement a functor method `(est::Subtype)(point::AbstractVector{<:Real}) -> T`.
"""
abstract type AbstractEstimator{T} end

"""
    (est::AbstractEstimator{T})(points::AbstractVector{<:AbstractVector})

Interpolate at a list of points using multithreading.
"""
function (est::AbstractEstimator{T})(points::AbstractVector{<:AbstractVector}) where T
    results = Vector{T}(undef, length(points))
    Threads.@threads for i in eachindex(points)
        results[i] = est(points[i])
    end
    return results
end


"""
    (est::AbstractEstimator{T})(grid::Tuple)

Interpolate on a 3D grid defined by (xs, ys, zs).
"""
function (est::AbstractEstimator{T})(grid::Tuple) where T
    xs, ys, zs = grid

    results = Array{T,3}(undef, length(xs), length(ys), length(zs))

    Threads.@threads for i in eachindex(xs)
        x = xs[i]
        @inbounds for (j, y) in pairs(ys)
            @inbounds for (k, z) in pairs(zs)
                results[i, j, k] = est(SVector(x, y, z))
            end
        end
    end
    return results
end


# -----------------------------------------------------------------------------
# DensityEstimator
# -----------------------------------------------------------------------------

"""
    DensityEstimator

A struct containing the necessary components for density estimation.
Output type: `Float64`
"""
struct DensityEstimator <: AbstractEstimator{Float64}
    bvh::BoundingVolumeHierarchy
    triangulation::Triangulation3D
    tetrahedra::Matrix{Int}
end

"""
    DensityEstimator(points, weights::Vector, depth=9)
    DensityEstimator(points, depth=9)

Build a DTFE estimator from a point cloud.
"""
function DensityEstimator(points, weights::Vector, depth::Int=9)
    coords, tets = tessellate(points)
    triangulation = Triangulation3D(points, tets, weights)  # tets is M×4, rows are tets

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    return DensityEstimator(bvh, triangulation, tets)
end

function DensityEstimator(points, depth::Int=9)
    coords, tets = tessellate(points)
    triangulation = Triangulation3D(points, tets)  # tets is M×4, rows are tets

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    return DensityEstimator(bvh, triangulation, tets)
end

"""
    (est::DensityEstimator)(point::AbstractVector{<:Real})

Interpolate density at a single point using barycentric coordinates.
"""
function (est::DensityEstimator)(point::AbstractVector{<:Real})
    dtfe(point, est.bvh, est.tetrahedra, est.triangulation)
end


# -----------------------------------------------------------------------------
# VelocityEstimator
# -----------------------------------------------------------------------------

"""
    VelocityEstimator

A struct containing the necessary components for velocity field estimation.
Output type: `SVector{3, Float64}`
"""
struct VelocityEstimator <: AbstractEstimator{SVector{3,Float64}}
    bvh::BoundingVolumeHierarchy
    triangulation::Triangulation3D
    tetrahedra::Matrix{Int}
    velocities::Vector{SVector{3,Float64}}
end

"""
    VelocityEstimator(points, velocities::Vector, depth=9)

Build a DTFE velocity estimator from a point cloud and associated velocities.
"""
function VelocityEstimator(points, velocities::Vector, depth::Int=9)
    length(velocities) == length(points) || throw(ArgumentError("length of velocities must match length of points"))
    coords, tets = tessellate(points)
    triangulation = Triangulation3D(points, tets) # Reusing existing Triangulation3D for topology

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    # Convert velocities to SVector if they aren't already
    vels = [SVector{3,Float64}(v) for v in velocities]

    return VelocityEstimator(bvh, triangulation, tets, vels)
end


"""
    VelocityEstimator(est::DensityEstimator, velocities::Vector)

Build a VelocityEstimator reusing the topology from a DensityEstimator.
"""
function VelocityEstimator(est::DensityEstimator, velocities::Vector)
    length(velocities) == length(est.triangulation.points) || throw(ArgumentError("length of velocities must match length of points"))
    vels = [SVector{3,Float64}(v) for v in velocities]
    return VelocityEstimator(est.bvh, est.triangulation, est.tetrahedra, vels)
end

"""
    (est::VelocityEstimator)(point::AbstractVector{<:Real})

Interpolate velocity at a single point.
"""
function (est::VelocityEstimator)(point::AbstractVector{<:Real})
    coords = est.triangulation.points

    i = findId(point, coords, est.tetrahedra, est.bvh)
    i === nothing && return SVector(0.0, 0.0, 0.0)

    @inbounds idx1 = est.tetrahedra[i, 1]
    @inbounds idx2 = est.tetrahedra[i, 2]
    @inbounds idx3 = est.tetrahedra[i, 3]
    @inbounds idx4 = est.tetrahedra[i, 4]

    @inbounds v1 = coords[idx1]
    @inbounds v2 = coords[idx2]
    @inbounds v3 = coords[idx3]
    @inbounds v4 = coords[idx4]

    @inbounds vel1 = est.velocities[idx1]
    @inbounds vel2 = est.velocities[idx2]
    @inbounds vel3 = est.velocities[idx3]
    @inbounds vel4 = est.velocities[idx4]

    dX = SMatrix{3,3}(hcat(v2 - v1, v3 - v1, v4 - v1))
    invM = inv(dX)

    diff = point - v1
    λ234 = invM * diff
    λ1 = 1.0 - sum(λ234)

    return λ1 * vel1 + λ234[1] * vel2 + λ234[2] * vel3 + λ234[3] * vel4
end

"""
    velocity(est::VelocityEstimator, point::AbstractVector{<:Real})

Interpolate velocity at a single point (alias for functor).
"""
function velocity(est::VelocityEstimator, point::AbstractVector{<:Real})
    est(point)
end

"""
    velocityGradient(est::VelocityEstimator, point::AbstractVector{<:Real})

Calculate interpolated velocity and its gradients (divergence, shear, vorticity) at a point.

# Returns
A tuple `(velocity, divergence, shear, vorticity)` where:
- `velocity` - Interpolated velocity vector
- `divergence` (∇⋅v)
- `shear` - Symmetric traceless shear tensor
- `vorticity` (∇×v)

Returns zero values if the point is outside the tessellation.
"""
function velocityGradient(est::VelocityEstimator, point::AbstractVector{<:Real})
    coords = est.triangulation.points

    i = findId(point, coords, est.tetrahedra, est.bvh)
    i === nothing && return (SVector(0.0, 0.0, 0.0), 0.0, zero(SMatrix{3,3,Float64}), SVector(0.0, 0.0, 0.0))

    @inbounds idx1 = est.tetrahedra[i, 1]
    @inbounds idx2 = est.tetrahedra[i, 2]
    @inbounds idx3 = est.tetrahedra[i, 3]
    @inbounds idx4 = est.tetrahedra[i, 4]

    @inbounds v1 = coords[idx1]
    @inbounds v2 = coords[idx2]
    @inbounds v3 = coords[idx3]
    @inbounds v4 = coords[idx4]

    @inbounds vel1 = est.velocities[idx1]
    @inbounds vel2 = est.velocities[idx2]
    @inbounds vel3 = est.velocities[idx3]
    @inbounds vel4 = est.velocities[idx4]


    # ΔX matrix columns: x2-x1, x3-x1, x4-x1
    dX = SMatrix{3,3}(hcat(v2 - v1, v3 - v1, v4 - v1))
    invM = inv(dX)

    diff = point - v1
    λ234 = invM * diff
    λ1 = 1.0 - sum(λ234)

    # v(x) = Σ λ_i v_i
    vInterp = λ1 * vel1 + λ234[1] * vel2 + λ234[2] * vel3 + λ234[3] * vel4

    dV = SMatrix{3,3}(hcat(vel2 - vel1, vel3 - vel1, vel4 - vel1))

    J = dV * invM # Velocity Gradient Tensor (∇v)

    div = tr(J)

    # Shear σ = 0.5(J + J^T) - (1/3)θI
    S = 0.5 * (J + transpose(J))
    shear = S - (div / 3.0) * I

    # Vorticity ω = ∇×v
    vorticity = SVector(
        J[3, 2] - J[2, 3],
        J[1, 3] - J[3, 1],
        J[2, 1] - J[1, 2]
    )

    return (vInterp, div, shear, vorticity)
end

# -----------------------------------------------------------------------------
# Shared Kernel
# -----------------------------------------------------------------------------

"""
    dtfe(point, bvh, tetrahedra, tessellation)

Interpolate density at a single point using barycentric coordinates.

# Returns
- Density value at point, or 0 if point is outside tessellation.
"""
@views function dtfe(point, bvh, tetrahedra, tessellation)
    coords = tessellation.points

    i = findId(point, coords, tetrahedra, bvh) # passes arrays without allocating to speed up
    i === nothing && return 0.0


    @inbounds idx1 = tetrahedra[i, 1]
    @inbounds idx2 = tetrahedra[i, 2]
    @inbounds idx3 = tetrahedra[i, 3]
    @inbounds idx4 = tetrahedra[i, 4]

    @inbounds v1 = coords[idx1]
    @inbounds v2 = coords[idx2]
    @inbounds v3 = coords[idx3]
    @inbounds v4 = coords[idx4]

    @inbounds rho1 = tessellation.rhoStar[idx1]
    @inbounds rho2 = tessellation.rhoStar[idx2]
    @inbounds rho3 = tessellation.rhoStar[idx3]
    @inbounds rho4 = tessellation.rhoStar[idx4]

    dX = SMatrix{3,3,Float64}(hcat(v2 - v1, v3 - v1, v4 - v1))
    invM = inv(dX)

    diff = point - v1
    λ234 = invM * diff
    λ1 = 1.0 - sum(λ234)

    interpolation = λ1 * rho1 + λ234[1] * rho2 + λ234[2] * rho3 + λ234[3] * rho4

    return interpolation
end


