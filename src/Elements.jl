# Elements - Core data structures for DTFE

"""3D point type using StaticArrays for performance."""
const Point3 = SVector{3,Float64}



"""
    Triangulation3D

Result of Delaunay tessellation with per-vertex density estimates.
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

function Triangulation3D(points::Vector{Point3}, tets::AbstractMatrix, weights::Vector)
    length(weights) == length(points) || throw(ArgumentError("length of weights must match length of points"))

    # DTFE density formula: ρ*ᵢ = (d+1) * wᵢ / Σⱼ Vⱼ

    nPoints = length(points)
    rhos = zeros(nPoints)

    #accumulate contiguous volume for each vertex
    @inbounds for tet in eachrow(tets)
        pos = points[tet]
        vol = computeVolume(pos)
        for i in tet
            rhos[i] += vol
        end
    end

    #compute density using DTFE formula in-place
    @inbounds for i in 1:nPoints
        if rhos[i] > 0
            rhos[i] = 4.0 * weights[i] / rhos[i]
        end
    end

    return Triangulation3D(points, rhos)
end
