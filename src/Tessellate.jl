# Tessellate - Delaunay tetrahedralization via TetGen
#
# Produces 3D tessellation from point clouds.

"""Convert vector of Point3 to matrix for TetGen."""
pointsToMatrix(points::Vector{Point3}) = reduce(hcat, points)'


"""
    tessellate(points::Vector{Point3})

Tetrahedralize a 3D point cloud.

# Returns
- `coords::Matrix{Float64}`: Coordinate matrix [3×N]
- `tetrahedra::Matrix{Int}`: Tetrahedron indices [M×4]
"""
function tessellate(points::Vector{Point3})
    coords3 = pointsToMatrix(points)

    meshdata = TetGen.RawTetGenIO{Float64}()
    meshdata.pointlist = coords3'

    result = TetGen.tetrahedralize(meshdata, "Q")

    tets = result.tetrahedronlist'
    coords = result.pointlist
    return coords, tets
end

function tessellate(points::Matrix)
    meshdata = TetGen.RawTetGenIO{Float64}()
    meshdata.pointlist = points

    result = TetGen.tetrahedralize(meshdata, "Q")

    tets = result.tetrahedronlist'
    coords = result.pointlist
    return coords, tets
end


