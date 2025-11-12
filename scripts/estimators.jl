module Estimators

using ..Elements, ..Tesselate, ..Bvh, ..Searchers
using StaticArrays
using LinearAlgebra
using KernelAbstractions
using CUDA

export standardEstimator, DTFE



tesselate = Tesselate.tesselate # this line should be unnecessary, but it doesn't work without it

function standardEstimator(points, weights, depth = 9)

    coords, tets = tesselate(points)
    tes = Triangulation3D(points,tets',weights)

    simplices = coords[:,tets]
    bvh = BVH(simplices,depth)

    return bvh,tes,tets
end


function standardEstimator(points,depth = 9)

    coords, tets = tesselate(points)
    tes = Triangulation3D(points,tets')

    simplices = coords[:,tets]
    bvh = BVH(simplices,depth)

    return bvh,tes,tets
end


function DTFE(point,bvh,tetrahedra,tesselation)

    coords = tesselation.points
    simplices = coords[tetrahedra]

    i = findID(point,simplices,bvh)
    if i == nothing
        return 0
    end

    tet = tetrahedra[i,:]
    simp = coords[tet]
    rhos = tesselation.ρStar[tet]

    delRho=invertClassic(rhos,simp)

    interpolation = rhos[1] + dot((point - simp[1]),delRho)

    return interpolation
end

@inline function invertClassic(rhos,simplex) 
    r = rhos[2:end] .- rhos[1]

    v1, v2, v3, v4 = simplex[1], simplex[2], simplex[3], simplex[4]

    a = v2 - v1 
    b = v3 - v1 
    c = v4 - v1

    mat = SMatrix{3,3}([a;b;c])
    
    return inv(mat)*r
end

@inline function invertGPU(rhos, simplex)
    r = rhos[2:end] .- rhos[1]
    a = simplex[2] - simplex[1]
    b = simplex[3] - simplex[1]
    c = simplex[4] - simplex[1]
    mat = SMatrix{3,3}([a'; b'; c'])
    return inv(mat) * r
end

@kernel function interpolateGPU!(points, tets, coords, rhoStars, out)
    i = @index(Global)
    if i <= size(points, 1)
        tet  = tets[i, :]
        simp = coords[tet, :]
        rhos = rhoStars[tet]
        delRho = invertGPU(rhos, simp)
        out[i] = rhos[1] + dot(points[i, :] .- simp[1, :], delRho)
    end
end

function DTFE(points::Matrix, bvh, tetrahedra, tesselation)

    ids = [findID(p, tesselation.points[tetrahedra], bvh) for p in eachrow(points)]

    valid = findall(!isnothing, ids)
    cleanPoints = points[:,valid]
    cleanIds = getindex.(ids[valid]) 

    coords = tesselation.points
    rhoStar = tesselation.ρStar
    tets = tetrahedra[cleanIds, :]

    cuPoints = CuArray(SVector{3,Float64}(cleanPoints))
    cuTets = CuArray(tets)
    cuCoords = CuArray(coords)
    cuRhoStar = CuArray(rhoStar)
    cuOut = CuArray(zeros(Float32, length(cleanPoints)))

    n = length(cleanPoints)
    interpolateGPU!(CUDADevice())(cuPoints, cuTets, cuCoords, cuRhoStar, cuOut; ndrange=n)

    results = Array(cuOut)

    out = zeros(Float32, length(points))
    out[valid] = results

    return out
end


end