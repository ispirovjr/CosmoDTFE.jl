@inline function recursiveSearch(point, bvh::BoundingVolumeHierarchy)
    return recursiveSearch(Point3(point), bvh.tree, bvh.bbox)
end

@inline recursiveSearch(::Point3, leaf::BVHLeaf, ::SearchBox3D) = leaf.data

@inline function recursiveSearch(
    point::Point3,
    tree::BVHNode,
    searchBox::SearchBox3D,
)
    ax = tree.depth % 3 + 1
    line = (searchBox[ax, 2] + searchBox[ax, 1]) / 2

    if point[ax] < line
        newBox = setindex(searchBox, line, ax, 2)
        return recursiveSearch(point, tree.leftChild, newBox)
    end

    newBox = setindex(searchBox, line, ax, 1)
    return recursiveSearch(point, tree.rightChild, newBox)
end

"""
    findSimplex(point, simplices, bvh)

Find the simplex containing `point`. Returns `nothing` when no simplex contains
the point. 
"""
function findSimplex(point, simplices::AbstractMatrix, bvh::BoundingVolumeHierarchy)
    indices = recursiveSearch(point, bvh)
    simplexNeighborhood = simplices[:, indices]

    idx = earlyStopSearch(point, simplexNeighborhood)
    idx === nothing && return nothing

    return simplexNeighborhood[:, idx]
end

"""
    findId(point, simplices, bvh)

Find the first simplex index containing `point`. 
"""
function findId(point::AbstractVector, simplices::AbstractMatrix, bvh::BoundingVolumeHierarchy)
    p = Point3(point)
    return findId(p, simplices, bvh)
end

function findId(point::SVector{3,Float64}, simplices::AbstractMatrix, bvh::BoundingVolumeHierarchy)
    indices = recursiveSearch(point, bvh)
    simplexNeighborhood = @view simplices[:, indices]
    idx = earlyStopSearch(point, simplexNeighborhood)

    idx === nothing && return nothing
    return indices[idx]
end

"""
    findId(point, coords, tetrahedra, bvh)

Memory-efficient point lookup that avoids materializing the full simplex list.
Returns the first containing tetrahedron id, or `nothing`.
"""
function findId(point::AbstractVector, coords::AbstractVector, tetrahedra::AbstractMatrix{<:Integer}, bvh::BoundingVolumeHierarchy)
    p = Point3(point)
    indices = recursiveSearch(p, bvh)

    @inbounds for tetId in indices
        v1 = coords[tetrahedra[1, tetId]]
        v2 = coords[tetrahedra[2, tetId]]
        v3 = coords[tetrahedra[3, tetId]]
        v4 = coords[tetrahedra[4, tetId]]

        if intersection3D(p, (v1, v2, v3, v4))
            return tetId
        end
    end
    return nothing
end

"""
    findAllIds(point, coords, tetrahedra, bvh)

Return all tetrahedra in the BVH leaf neighborhood that contain `point`.
This is intended for warped phase-space tessellations where tetrahedra can
overlap in evaluation coordinates.
"""
function findAllIds(point::AbstractVector, coords::AbstractVector, tetrahedra::AbstractMatrix{<:Integer}, bvh::BoundingVolumeHierarchy)
    staticPoint = Point3(point)
    indices = recursiveSearch(staticPoint, bvh)
    matches = Int[]

    @inbounds for tetrahedronId in indices
        vertexOne = coords[tetrahedra[1, tetrahedronId]]
        vertexTwo = coords[tetrahedra[2, tetrahedronId]]
        vertexThree = coords[tetrahedra[3, tetrahedronId]]
        vertexFour = coords[tetrahedra[4, tetrahedronId]]

        if intersection3D(staticPoint, (vertexOne, vertexTwo, vertexThree, vertexFour))
            push!(matches, tetrahedronId)
        end
    end

    return matches
end

function earlyStopSearch(point::SVector{3,Float64}, simplices::AbstractMatrix)
    for (i, simplex) in enumerate(eachcol(simplices))
        if intersection3D(point, simplex)
            return i
        end
    end
    return nothing
end

function earlyStopSearch(point::AbstractVector, simplices::AbstractMatrix)
    return earlyStopSearch(Point3(point), simplices)
end

@inline function intersection3D(point::SVector{3,Float64}, simplex::AbstractMatrix{<:Real})
    @inbounds begin
        v1 = Point3(simplex[1, 1], simplex[2, 1], simplex[3, 1])
        v2 = Point3(simplex[1, 2], simplex[2, 2], simplex[3, 2])
        v3 = Point3(simplex[1, 3], simplex[2, 3], simplex[3, 3])
        v4 = Point3(simplex[1, 4], simplex[2, 4], simplex[3, 4])
    end
    return intersection3D(point, (v1, v2, v3, v4))
end

@inline function intersection3D(point::SVector{3,Float64}, simplex)
    @inbounds begin
        v1, v2, v3, v4 = simplex[1], simplex[2], simplex[3], simplex[4]
        a = v2 - v1
        b = v3 - v1
        c = v4 - v1
        r = point - v1

        detA = dot(a, cross(b, c))
        lambda2 = dot(r, cross(b, c)) / detA
        lambda3 = dot(r, cross(c, a)) / detA
        lambda4 = dot(r, cross(a, b)) / detA
    end

    s = lambda2 + lambda3 + lambda4
    return (lambda2 >= 0) & (lambda3 >= 0) & (lambda4 >= 0) & (s <= 1)
end

@inline function intersection3D(point::AbstractVector, simplex)
    return intersection3D(Point3(point), simplex)
end
