module Estimators

using ..Elements, ..Tesselate, ..Bvh, ..Searchers
using StaticArrays
using LinearAlgebra
using KernelAbstractions
using CUDA

export standardEstimator, DTFE, DTFEMultiThread



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


@views function DTFE(point,bvh,tetrahedra,tesselation)

    coords = tesselation.points
    simplices = coords[tetrahedra]

    i = findID(point,simplices,bvh)
    i === nothing && return 0

    tet = tetrahedra[i,:]
    simp = coords[tet]
    rhos = tesselation.ρStar[tet]

#    delRho=invertClassic(rhos,simp)

    v1 = simp[1]
    invM = inv(SMatrix{3,3}(hcat(simp[2]-v1, simp[3]-v1, simp[4]-v1)))
    diff = point - v1
    λ234 = invM * diff
    λ1   = 1 - sum(λ234)
    interpolation = λ1*rhos[1] + λ234[1]*rhos[2] + λ234[2]*rhos[3] + λ234[3]*rhos[4]

    #interpolation = rhos[1] + dot((point - simp[1]),delRho)

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

function DTFECuda(points::Matrix, bvh, tetrahedra, tesselation) #Abandon until we figure if we can use GPU in kapteyn

    ids = [findID(pt, tesselation.points[tetrahedra], bvh) for pt in points]

    valid = findall(!isnothing, ids)
    cleanPoints = points[valid]
    cleanIds = getindex.(ids[valid]) 

    coords = tesselation.points
    rhoStar = tesselation.ρStar
    tets = tetrahedra[cleanIds, :]

    cuPoints = [CuArray(SVector{3,Float64}(cleanPoint)) for cleanPoint in eachcol(cleanPoints)] 
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

# we multiprocess over the first dimension, so for heterogenious sampling, prioritize x size
function DTFEMultiThread(points, bvh, tetrahedra, tesselation) #TODO see if we can elegantly decide which dim to mutithread

    dens = zeros(size(points[1],1), size(points[2],1), size(points[3],1))


    xs,ys,zs = points

    Threads.@threads for i in eachindex(xs)
    x = xs[i]
    @inbounds for (j,y) in pairs(ys) 
            @inbounds for (k,z) in pairs(zs)
                dens[i,j,k] = DTFE([x,y,z], bvh, tetrahedra, tesselation) 
            end
        end
    end
    return dens
end




end