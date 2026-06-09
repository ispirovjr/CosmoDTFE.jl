"""
    AbstractEstimator{T}

Abstract base type for DTFE estimators that return values of type `T`.
Subtypes implement `(est::Subtype)(point::AbstractVector{<:Real})`.
"""
abstract type AbstractEstimator{T} end

"""
    (est::AbstractEstimator{T})(points)

Evaluate an estimator at a vector of points using Julia threads.
"""
function (est::AbstractEstimator{T})(points::AbstractVector{<:AbstractVector}) where T
    results = Vector{T}(undef, length(points))
    Threads.@threads for i in eachindex(points)
        results[i] = est(points[i])
    end
    return results
end

"""
    (est::AbstractEstimator{T})(grid)

Evaluate an estimator on a 3D grid represented by `(xs, ys, zs)`.
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

@inline function barycentricWeights(point::AbstractVector{<:Real}, v1, v2, v3, v4)
    p = Point3(point)
    dX = SMatrix{3,3,Float64}(hcat(v2 - v1, v3 - v1, v4 - v1))
    lambda234 = inv(dX) * (p - v1)
    lambda1 = 1.0 - sum(lambda234)
    return SVector(lambda1, lambda234[1], lambda234[2], lambda234[3])
end

@inline function interpolateSimplex(lambda, val1, val2, val3, val4)
    return lambda[1] * val1 + lambda[2] * val2 + lambda[3] * val3 + lambda[4] * val4
end

@inline function simplexData(data, tets, tetId::Integer)
    @inbounds return (
        data[tets[tetId, 1]],
        data[tets[tetId, 2]],
        data[tets[tetId, 3]],
        data[tets[tetId, 4]],
    )
end

function buildEstimatorGeometry(points, depth::Int)
    ps = toPoint3Vector(points)
    coords, tets = tessellate(ps)
    tessPoints = matrixColumnsToPoints(coords)

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)
    return tessPoints, tets, bvh
end

"""
    DensityEstimator

DTFE density estimator. Calling the estimator at a point returns a `Float64`
density value, or `0.0` outside the tessellated volume.
"""
struct DensityEstimator <: AbstractEstimator{Float64}
    bvh::BoundingVolumeHierarchy
    triangulation::Triangulation3D
    tetrahedra::Matrix{Int}
end

"""
    DensityEstimator(points, weights; depth=9)
    DensityEstimator(points, weights, depth)
    DensityEstimator(points; depth=9)
    DensityEstimator(points, depth)

Build a DTFE density estimator from point positions and optional per-point
weights. Unweighted density uses unit weights.
"""
function DensityEstimator(points, weights::AbstractVector; depth::Int=9)
    pts, tets, bvh = buildEstimatorGeometry(points, depth)
    triangulation = Triangulation3D(pts, tets, weights)
    return DensityEstimator(bvh, triangulation, tets)
end

DensityEstimator(points, weights::AbstractVector, depth::Int) = DensityEstimator(points, weights; depth=depth)

function DensityEstimator(points; depth::Int=9)
    pointCount = length(points)
    weights = ones(Float64, pointCount)
    return DensityEstimator(points, weights; depth=depth)
end

DensityEstimator(points, depth::Int) = DensityEstimator(points; depth=depth)

"""
    (est::DensityEstimator)(point)

Interpolate density at a single point.
"""
function (est::DensityEstimator)(point::AbstractVector{<:Real})
    return dtfe(point, est.bvh, est.tetrahedra, est.triangulation)
end

"""
    VelocityEstimator

DTFE velocity estimator. Calling the estimator at a point returns an
`SVector{3,Float64}` velocity, or a zero vector outside the tessellated volume.
"""
struct VelocityEstimator <: AbstractEstimator{SVector{3,Float64}}
    bvh::BoundingVolumeHierarchy
    triangulation::Triangulation3D
    tetrahedra::Matrix{Int}
    velocities::Vector{SVector{3,Float64}}
end

"""
    VelocityEstimator(points, velocities; depth=9)
    VelocityEstimator(points, velocities, depth)

Build a DTFE velocity estimator from point positions and per-point velocities.
"""
function VelocityEstimator(points, velocities::AbstractVector; depth::Int=9)
    length(velocities) == length(points) || throw(ArgumentError("length of velocities must match length of points"))

    pts, tets, bvh = buildEstimatorGeometry(points, depth)
    triangulation = Triangulation3D(pts, tets)
    vels = [Point3(velocity) for velocity in velocities]

    return VelocityEstimator(bvh, triangulation, tets, vels)
end

VelocityEstimator(points, velocities::AbstractVector, depth::Int) = VelocityEstimator(points, velocities; depth=depth)

"""
    VelocityEstimator(est::DensityEstimator, velocities)

Build a velocity estimator reusing the tessellation from an existing density
estimator.
"""
function VelocityEstimator(est::DensityEstimator, velocities::AbstractVector)
    length(velocities) == length(est.triangulation.points) ||
        throw(ArgumentError("length of velocities must match length of points"))
    @warn "Generating velocity field from density field. If the density field was mass weighted, this estimates momenta, not velocities."
    vels = [Point3(velocity) for velocity in velocities]
    return VelocityEstimator(est.bvh, est.triangulation, est.tetrahedra, vels)
end

"""
    (est::VelocityEstimator)(point)

Interpolate velocity at a single point.
"""
function (est::VelocityEstimator)(point::AbstractVector{<:Real})
    tetId = findId(point, est.triangulation.points, est.tetrahedra, est.bvh)
    tetId === nothing && return Point3(0.0, 0.0, 0.0)

    verts = simplexData(est.triangulation.points, est.tetrahedra, tetId)
    vels = simplexData(est.velocities, est.tetrahedra, tetId)
    lambda = barycentricWeights(point, verts...)

    return interpolateSimplex(lambda, vels...)
end

"""
    velocity(est::VelocityEstimator, point)

Interpolate velocity at a single point.
"""
velocity(est::VelocityEstimator, point::AbstractVector{<:Real}) = est(point)

"""
    velocityGradient(est::VelocityEstimator, point)

Calculate interpolated velocity and its local gradient decomposition.

Returns `(velocity, divergence, shear, vorticity)`. All values are zero when
the point is outside the tessellated volume.
"""
function velocityGradient(est::VelocityEstimator, point::AbstractVector{<:Real})
    tetId = findId(point, est.triangulation.points, est.tetrahedra, est.bvh)
    tetId === nothing && return (Point3(0.0, 0.0, 0.0), 0.0, zero(SMatrix{3,3,Float64}), Point3(0.0, 0.0, 0.0))

    v1, v2, v3, v4 = simplexData(est.triangulation.points, est.tetrahedra, tetId)
    vel1, vel2, vel3, vel4 = simplexData(est.velocities, est.tetrahedra, tetId)

    dX = SMatrix{3,3,Float64}(hcat(v2 - v1, v3 - v1, v4 - v1))
    invM = inv(dX)

    lambda = barycentricWeights(point, v1, v2, v3, v4)
    vInterp = interpolateSimplex(lambda, vel1, vel2, vel3, vel4)

    dV = SMatrix{3,3,Float64}(hcat(vel2 - vel1, vel3 - vel1, vel4 - vel1))
    J = dV * invM

    div = tr(J)
    S = 0.5 * (J + transpose(J))
    I3 = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    shear = S - (div / 3.0) * I3

    vorticity = Point3(
        J[3, 2] - J[2, 3],
        J[1, 3] - J[3, 1],
        J[2, 1] - J[1, 2],
    )

    return vInterp, div, shear, vorticity
end

"""
    dtfe(point, bvh, tetrahedra, tessellation)

Interpolate density at a single point using barycentric coordinates. Returns
`0.0` when the point is outside the tessellated volume.
"""
function dtfe(point, bvh, tetrahedra, tessellation)
    tetId = findId(point, tessellation.points, tetrahedra, bvh)
    tetId === nothing && return 0.0

    verts = simplexData(tessellation.points, tetrahedra, tetId)
    rhos = simplexData(tessellation.rhoStar, tetrahedra, tetId)
    lambda = barycentricWeights(point, verts...)

    return interpolateSimplex(lambda, rhos...)
end
