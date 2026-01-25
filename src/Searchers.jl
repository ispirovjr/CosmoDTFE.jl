# Searchers - Spatial search algorithms using BVH
#
# Provides fast point-in-tetrahedron queries.

function recursiveSearch(point, leaf::BVHLeaf, bbox::Matrix)
    return leaf.data
end

function recursiveSearch(point, tree::BVHNode, bbox::Matrix)
    ax = tree.depth % 3 + 1

    newBox = copy(bbox)

    line = (bbox[ax, 2] + bbox[ax, 1]) / 2
    if point[ax] < line
        newBox[ax, 2] = line
        return recursiveSearch(point, tree.leftChild, newBox)
    end
    newBox[ax, 1] = line
    return recursiveSearch(point, tree.rightChild, newBox)
end


"""
    findSimplex(point, simplices, bvh)

Find the simplex containing the given point.
"""
function findSimplex(point, simplices, bvh::BoundingVolumeHierarchy)
    indices = recursiveSearch(point, bvh.tree, bvh.bbox)

    simplexNeighborhood = simplices[indices, :]

    idx = earlyStopSearch(point, simplexNeighborhood)

    return simplexNeighborhood[idx]
end

"""
    findId(point, simplices, bvh)

Find the index of the simplex containing the given point.
Returns `nothing` if point is outside all simplices.
"""
function findId(point, simplices, bvh::BoundingVolumeHierarchy)
    indices = recursiveSearch(point, bvh.tree, bvh.bbox)

    simplexNeighborhood = simplices[indices, :]

    idx = earlyStopSearch(point, simplexNeighborhood)

    if idx === nothing
        return nothing
    end

    return indices[idx]
end

function findId(point, simplices::Array, bvh::BoundingVolumeHierarchy)
    indices = recursiveSearch(point, bvh.tree, bvh.bbox)

    simplexNeighborhood = simplices[:, indices, :]

    idx = earlyStopSearch(point, simplexNeighborhood)

    if idx === nothing
        return nothing
    end

    return indices[idx]
end

function findId(point, simplices::Array{SVector{3,Float64},2}, bvh::BoundingVolumeHierarchy)
    indices = recursiveSearch(point, bvh.tree, bvh.bbox)

    simplexNeighborhood = simplices[indices, :]

    idx = earlyStopSearch(point, simplexNeighborhood)

    if idx === nothing
        return nothing
    end

    return indices[idx]
end


function earlyStopSearch(p::Vector, simplices::Matrix)
    for (i, s) in pairs(eachrow(simplices))
        if intersection3D(p, s)
            return i
        end
    end
end

# Optimal path: StaticArrays with Cramer's rule
# NOTE: Using bitwise & instead of && for performance (avoids branch prediction overhead)
@inline function intersection3D(p::SVector{3,Float64}, simplex::SMatrix{4,3,Float64})
    @inbounds begin
        v1, v2, v3, v4 = simplex[1, :], simplex[2, :], simplex[3, :], simplex[4, :]
        a = v2 - v1
        b = v3 - v1
        c = v4 - v1
        r = p - v1

        detA = dot(a, cross(b, c))
        x1 = dot(r, cross(b, c)) / detA
        x2 = dot(r, cross(c, a)) / detA
        x3 = dot(r, cross(a, b)) / detA
    end

    s = x1 + x2 + x3
    return (x1 >= 0) & (x2 >= 0) & (x3 >= 0) & (s <= 1)
end

# Float vector with generic simplex type
@inline function intersection3D(point::Vector{Float64}, simplex)
    @inbounds begin
        p = SVector{3,Float64}(point)

        v1, v2, v3, v4 = simplex[1], simplex[2], simplex[3], simplex[4]
        a = v2 - v1
        b = v3 - v1
        c = v4 - v1
        r = p - v1

        detA = dot(a, cross(b, c))
        x1 = dot(r, cross(b, c)) / detA
        x2 = dot(r, cross(c, a)) / detA
        x3 = dot(r, cross(a, b)) / detA
    end

    s = x1 + x2 + x3
    return (x1 >= 0) & (x2 >= 0) & (x3 >= 0) & (s <= 1)
end

# Matrix input: convert to StaticArrays
function intersection3D(p::Vector{Float64}, simplex::Matrix)
    sP = SVector{3}(p)
    sSimp = SMatrix{4,3}(simplex)
    return intersection3D(sP, sSimp)
end
