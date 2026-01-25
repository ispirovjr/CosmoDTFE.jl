# Elements - Core data structures for DTFE
#
# Provides Point3, Tetrahedron, and Triangulation3D types.

"""3D point type using StaticArrays for performance."""
const Point3 = SVector{3,Float64}


"""
    Tetrahedron

A tetrahedron defined by 4 vertices with precomputed volume.

# Fields
- `verts::NTuple{4, Point3}`: The 4 vertices
- `vol::Float64`: Precomputed volume
"""
struct Tetrahedron
    verts::NTuple{4,Point3}
    vol::Float64
end


"""
    Triangulation3D

Result of Delaunay tessellation with per-vertex density estimates.

# Fields
- `points::Vector{Point3}`: Original point cloud
- `rhoStar::Vector{Float64}`: Density estimate at each point
"""
struct Triangulation3D
    points::Vector{Point3}
    rhoStar::Vector{Float64}
end


"""
    computeVolume(verts)

Compute volume of tetrahedron from 4 vertices.
"""
function computeVolume(verts)
    v1, v2, v3, v4 = verts
    a = v2 - v1
    b = v3 - v1
    c = v4 - v1
    return abs(dot(a, cross(b, c))) / 6
end

function computeVolume(verts::Matrix)
    v1, v2, v3, v4 = eachcol(verts)
    a = v2 - v1
    b = v3 - v1
    c = v4 - v1
    return abs(dot(a, cross(b, c))) / 6
end


function Tetrahedron(verts::NTuple{4,Point3})
    vol = computeVolume(verts)
    return Tetrahedron(verts, vol)
end


function Triangulation3D(points::Vector{Point3}, tets::Matrix)
    weights = ones(length(points))
    return Triangulation3D(points, tets, weights)
end

function Triangulation3D(points::Vector{Point3}, tets::Matrix, weights::Vector)
    rhos = zeros(length(points))

    for tet in eachrow(tets)
        pos = points[tet]
        vol = computeVolume(pos)
        for i in tet
            rhos[i] += weights[i] / vol
        end
    end

    return Triangulation3D(points, rhos)
end
