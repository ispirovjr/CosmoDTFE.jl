module Elements

using LinearAlgebra
using StaticArrays

export Tetrahedron, Triangulation3D, point3,point2

const point3 = SVector{3,Float64}
const point2 = SVector{3,Float64}



struct Tetrahedron
    verts::NTuple{4, point3}    
    vol::Float64
end


struct Triangulation3D
    points::Vector{point3}
    density::Vector{Float64}
end


# --- Utility function to compute volume ---
function computeVolume(verts::NTuple{4, point3})
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


function Tetrahedron(verts::NTuple{4, point3})
    vol = computeVolume(verts)
    return Tetrahedron(verts, vol)
end


function Triangulation3D(points::Vector{point3},tets::Matrix)
    weights = ones(size(points),1)
    return Triangulation3D(points,tets,weights)
end

function Triangulation3D(points::Vector{point3},tets::Matrix,weights::Vector{Float64})
    rhos = zeros(size(points))
    
    for tet in eachrow(tets)
        pos = points[:]

    end
    
end

end