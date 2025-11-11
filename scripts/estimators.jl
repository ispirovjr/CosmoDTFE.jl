module Estimators

using ..Elements, ..Tesselate, ..Bvh, ..Searchers
using StaticArrays
using LinearAlgebra

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

function invertClassic(rhos,simplex) 
    r = rhos[2:end] .- rhos[1]

    v1, v2, v3, v4 = simplex[1], simplex[2], simplex[3], simplex[4]

    a = v2 - v1 
    b = v3 - v1 
    c = v4 - v1

    mat = SMatrix{3,3}([a;b;c])
    
    return inv(mat)*r
end

end