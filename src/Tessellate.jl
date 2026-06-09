"""
    pointsToColumns(points)

Convert a vector of `Point3` values into the `3 x nPoints` coordinate matrix
used by TetGen and the BVH code.
"""
pointsToColumns(points::Vector{Point3}) = reduce(hcat, points)

"""
    matrixColumnsToPoints(coords)

Convert a `3 x nPoints` coordinate matrix into a vector of `Point3` values.
"""
matrixColumnsToPoints(coords::AbstractMatrix) = [Point3(coords[:, pointId]) for pointId in axes(coords, 2)]

"""
    tessellate(points::Vector{Point3})

Tetrahedralize a 3D point cloud.

# Returns
- `coords::Matrix{Float64}`: coordinate matrix with size `3 x nPoints`
- `tetrahedra::Matrix{Int}`: tetrahedron indices with size `nTetrahedra x 4`
"""
function tessellate(points::Vector{Point3})
    meshData = TetGen.RawTetGenIO{Float64}()
    meshData.pointlist = pointsToColumns(points)

    result = TetGen.tetrahedralize(meshData, "Q")

    coords = result.pointlist
    tetrahedra = result.tetrahedronlist'
    return coords, tetrahedra
end

"""
    tessellate(points::Matrix)

Tetrahedralize a `3 x nPoints` coordinate matrix.
"""
function tessellate(points::Matrix)
    meshData = TetGen.RawTetGenIO{Float64}()
    meshData.pointlist = points

    result = TetGen.tetrahedralize(meshData, "Q")

    coords = result.pointlist
    tetrahedra = result.tetrahedronlist'
    return coords, tetrahedra
end
