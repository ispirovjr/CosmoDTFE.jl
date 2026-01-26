# Estimators - DTFE density estimation functions


"""
    standardEstimator(points, weights::Vector, depth=9)

Build a DTFE estimator from a point cloud with weights.

# Returns
- `bvh::BoundingVolumeHierarchy`: Spatial index for fast queries
- `triangulation::Triangulation3D`: Tessellation with density estimates
- `tetrahedra::Matrix{Int}`: Tetrahedron indices
"""
function standardEstimator(points, weights::Vector, depth::Int=9)
    coords, tets = tessellate(points)
    triangulation = Triangulation3D(points, tets, weights)  # tets is M×4, rows are tets

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    return bvh, triangulation, tets
end

"""
    standardEstimator(points, depth=9)

Build a DTFE estimator from a point cloud with uniform weights.
"""
function standardEstimator(points, depth::Int=9)
    coords, tets = tessellate(points)
    triangulation = Triangulation3D(points, tets)  # tets is M×4, rows are tets

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    return bvh, triangulation, tets
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


@inline function invertGpu(rhos, simplex)
    r = rhos[2:end] .- rhos[1]
    a = simplex[2] - simplex[1]
    b = simplex[3] - simplex[1]
    c = simplex[4] - simplex[1]
    mat = SMatrix{3,3}([a'; b'; c'])
    return inv(mat) * r
end

@kernel function interpolateGpu!(points, tets, coords, rhoStars, out)
    i = @index(Global)
    if i <= size(points, 1)
        tet = tets[i, :]
        simp = coords[tet, :]
        rhos = rhoStars[tet]
        delRho = invertGpu(rhos, simp)
        out[i] = rhos[1] + dot(points[i, :] .- simp[1, :], delRho)
    end
end

# GPU implementation (experimental - not fully tested)
function dtfeCuda(points::Matrix, bvh, tetrahedra, tessellation)
    ids = [findId(pt, tessellation.points[tetrahedra], bvh) for pt in points]

    valid = findall(!isnothing, ids)
    cleanPoints = points[valid]
    cleanIds = getindex.(ids[valid])

    coords = tessellation.points
    rhoStar = tessellation.rhoStar
    tets = tetrahedra[cleanIds, :]

    cuPoints = [CuArray(SVector{3,Float64}(cleanPoint)) for cleanPoint in eachcol(cleanPoints)]
    cuTets = CuArray(tets)
    cuCoords = CuArray(coords)
    cuRhoStar = CuArray(rhoStar)
    cuOut = CuArray(zeros(Float32, length(cleanPoints)))

    n = length(cleanPoints)
    interpolateGpu!(CUDADevice())(cuPoints, cuTets, cuCoords, cuRhoStar, cuOut; ndrange=n)

    results = Array(cuOut)

    out = zeros(Float32, length(points))
    out[valid] = results

    return out
end

"""
    dtfeMultiThread(points, bvh, tetrahedra, tessellation)

Multi-threaded field estimation on a 3D grid.

# Arguments
- `points`: Tuple of (xs, ys, zs) ranges for grid coordinates

# Returns
- `density::Array{Float64, 3}`: 3D density field
"""
function dtfeMultiThread(points, bvh, tetrahedra, tessellation)
    dens = zeros(length(points[1]), length(points[2]), length(points[3]))

    xs, ys, zs = points

    Threads.@threads for i in eachindex(xs)
        x = xs[i]
        @inbounds for (j, y) in pairs(ys)
            @inbounds for (k, z) in pairs(zs)
                dens[i, j, k] = dtfe([x, y, z], bvh, tetrahedra, tessellation)
            end
        end
    end
    return dens
end
