# Estimators - DTFE density estimation functions


"""
    DensityEstimator

A struct containing the necessary components for density estimation.
"""
struct DensityEstimator
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

"""
    (est::DensityEstimator)(points::AbstractVector{<:AbstractVector})

Interpolate density at a list of points using multithreading.
"""
function (est::DensityEstimator)(points::AbstractVector{<:AbstractVector})
    dens = zeros(length(points))
    Threads.@threads for i in eachindex(points)
        dens[i] = est(points[i])
    end
    return dens
end

"""
    (est::DensityEstimator)(grid::Tuple)

Interpolate density on a 3D grid defined by (xs, ys, zs).
"""
function (est::DensityEstimator)(grid::Tuple)
    xs, ys, zs = grid
    dens = zeros(length(xs), length(ys), length(zs))

    Threads.@threads for i in eachindex(xs)
        x = xs[i]
        @inbounds for (j, y) in pairs(ys)
            @inbounds for (k, z) in pairs(zs)
                # Use SVector for performance and correct dispatch
                dens[i, j, k] = est(SVector(x, y, z))
            end
        end
    end
    return dens
end


"""
    dtfe(point, bvh, tetrahedra, tessellation)

Interpolate density at a single point using barycentric coordinates.

# Returns
- Density value at point, or 0 if point is outside tessellation.
"""
@views function dtfe(point, bvh, tetrahedra, tessellation)
    coords = tessellation.points
    simplices = coords[tetrahedra]

    i = findId(point, simplices, bvh)
    i === nothing && return 0.0

    tet = tetrahedra[i, :]
    simp = coords[tet]
    rhos = tessellation.rhoStar[tet]

    v1 = simp[1]
    invM = inv(SMatrix{3,3}(hcat(simp[2] - v1, simp[3] - v1, simp[4] - v1)))
    diff = point - v1
    λ234 = invM * diff
    λ1 = 1 - sum(λ234)
    interpolation = λ1 * rhos[1] + λ234[1] * rhos[2] + λ234[2] * rhos[3] + λ234[3] * rhos[4]

    return interpolation
end

@inline function invertClassic(rhos, simplex)
    r = rhos[2:end] .- rhos[1]

    v1, v2, v3, v4 = simplex[1], simplex[2], simplex[3], simplex[4]

    a = v2 - v1
    b = v3 - v1
    c = v4 - v1

    mat = SMatrix{3,3}([a; b; c])

    return inv(mat) * r
end

# Deprecated/Legacy interface if needed
function dtfeMultiThread(points, bvh, tetrahedra, tessellation)
    est = DensityEstimator(bvh, tessellation, tetrahedra)
    return est(points)
end
