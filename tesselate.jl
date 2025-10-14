module Tesselate

using TetGen
using ..Elements

pointsToMatrix(points::Vector{point3}) = reduce(hcat, points)'


"""
tesselate3d(points::Matrix{Float64})
Tetrahedralizes a 3D point cloud.
Returns (coords, tetrahedra)
"""
function tesselate(points::Vector{point3})

    coords3 = pointsToMatrix(points)

    meshdata = TetGen.RawTetGenIO{Float64}()
    meshdata.pointlist = coords3'

    result = TetGen.tetrahedralize(meshdata, "Q")

    tets = result.tetrahedronlist'
    coords = result.pointlist
    return coords, tets
end

function tesselate(points::Matrix{Float64})
    meshdata = TetGen.RawTetGenIO()
    meshdata.pointlist = points'

    result = TetGen.tetrahedralize(meshdata, "Q")

    tets = result.tetrahedronlist'
    coords = result.pointlist
    return coords, tets
end

end