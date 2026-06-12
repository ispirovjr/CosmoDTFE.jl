"""3D point type using StaticArrays for performance."""
const Point3 = SVector{3,Float64}

"""Convert point-like inputs into a fresh vector of `Point3` values."""
toPoint3Vector(points::AbstractVector) = [Point3(point) for point in points]

"""
    Triangulation3D

Result of Delaunay tessellation with per-vertex DTFE density estimates.
"""
struct Triangulation3D{T}
    points::Vector{Point3}
    rhoStar::Vector{T}
end

"""
    computeVolume(verts)

Compute the volume of a tetrahedron from four vertices.
"""
function computeVolume(verts)
    v1, v2, v3, v4 = verts
    a = v2 - v1
    b = v3 - v1
    c = v4 - v1
    return abs(dot(a, cross(b, c))) / 6
end

function computeVolume(verts::AbstractMatrix)
    v1, v2, v3, v4 = eachcol(verts)
    a = v2 - v1
    b = v3 - v1
    c = v4 - v1
    return abs(dot(a, cross(b, c))) / 6
end

function Triangulation3D(points::Vector{Point3}, tets::AbstractMatrix)
    weights = ones(length(points))
    return Triangulation3D(points, tets, weights)
end

function Triangulation3D(points::Vector{Point3}, tets::AbstractMatrix, weights::AbstractVector)
    length(weights) == length(points) || throw(ArgumentError("length of weights must match length of points"))

    nPoints = length(points)
    rhoStar = zeros(nPoints)

    @inbounds for tet in eachrow(tets)
        pos = points[tet]
        vol = computeVolume(pos)
        for i in tet
            rhoStar[i] += vol
        end
    end

    @inbounds for i in 1:nPoints
        rhoStar[i] = 4.0 * weights[i] / rhoStar[i]
    end

    return Triangulation3D{Float64}(points, rhoStar)
end
