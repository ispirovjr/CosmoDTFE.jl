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


# ── Quickhull-based tessellation ──────────────────────────────────────────────

"""Convert Quickhull DelaunayHull facets to an M×4 Int matrix."""
function facetsToMatrix(hull)
    fs = collect(Quickhull.facets(hull))
    M = length(fs)
    tets = Matrix{Int}(undef, M, 4)
    @inbounds for (i, f) in enumerate(fs)
        tets[i, 1] = Int(f[1])
        tets[i, 2] = Int(f[2])
        tets[i, 3] = Int(f[3])
        tets[i, 4] = Int(f[4])
    end
    return tets
end


"""
    tessellateQH(points::Vector{Point3})

Tetrahedralize a 3D point cloud using Quickhull.

# Returns
- `coords::Matrix{Float64}`: Coordinate matrix [3×N]
- `tetrahedra::Matrix{Int}`: Tetrahedron indices [M×4]
"""
function tessellateQH(points::Vector{Point3})
    coords3 = pointsToMatrix(points)   # N×3
    coords  = coords3'                 # 3×N
    hull    = Quickhull.delaunay(coords)
    tets    = facetsToMatrix(hull)
    return coords, tets
end

function tessellateQH(points::Matrix)
    hull = Quickhull.delaunay(points)
    tets = facetsToMatrix(hull)
    return points, tets
end
