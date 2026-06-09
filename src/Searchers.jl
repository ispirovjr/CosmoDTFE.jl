function recursiveSearch(point, leaf::BVHLeaf, searchBox::Matrix)
    return leaf.data
end

function recursiveSearch(point, tree::BVHNode, searchBox::Matrix)
    ax = tree.depth % 3 + 1
    newBox = copy(searchBox)
    line = (searchBox[ax, 2] + searchBox[ax, 1]) / 2

    if point[ax] < line
        newBox[ax, 2] = line
        return recursiveSearch(point, tree.leftChild, newBox)
    end

    newBox[ax, 1] = line
    return recursiveSearch(point, tree.rightChild, newBox)
end

"""
    findSimplex(point, simplices, bvh)

Find the simplex containing `point`. Returns `nothing` when no simplex contains
the point.
"""
function findSimplex(point, simplices, bvh::BoundingVolumeHierarchy)
    indices = recursiveSearch(point, bvh.tree, bvh.bbox)
    simplexNeighborhood = simplices[indices, :]

    idx = earlyStopSearch(point, simplexNeighborhood)
    idx === nothing && return nothing

    return simplexNeighborhood[idx, :]
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
    indices = recursiveSearch(point, bvh.tree, bvh.bbox)
    simplexNeighborhood = @view simplices[indices, :]
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
    indices = recursiveSearch(p, bvh.tree, bvh.bbox)

    @inbounds for tetId in indices
        v1 = coords[tetrahedra[tetId, 1]]
        v2 = coords[tetrahedra[tetId, 2]]
        v3 = coords[tetrahedra[tetId, 3]]
        v4 = coords[tetrahedra[tetId, 4]]

        if intersection3D(p, (v1, v2, v3, v4))
            return tetId
        end
    end
    return nothing
end

function earlyStopSearch(point::SVector{3,Float64}, simplices::AbstractMatrix)
    for (i, simplex) in enumerate(eachrow(simplices))
        if intersection3D(point, simplex)
            return i
        end
    end
    return nothing
end

function earlyStopSearch(point::AbstractVector, simplices::AbstractMatrix)
    return earlyStopSearch(Point3(point), simplices)
end

@inline function intersection3D(point::SVector{3,Float64}, simplex::SMatrix{4,3,Float64})
    @inbounds begin
        v1 = Point3(simplex[1, 1], simplex[1, 2], simplex[1, 3])
        v2 = Point3(simplex[2, 1], simplex[2, 2], simplex[2, 3])
        v3 = Point3(simplex[3, 1], simplex[3, 2], simplex[3, 3])
        v4 = Point3(simplex[4, 1], simplex[4, 2], simplex[4, 3])
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

function intersection3D(point::AbstractVector, simplex::Matrix)
    p = Point3(point)
    s = SMatrix{4,3,Float64}(simplex)
    return intersection3D(p, s)
end
