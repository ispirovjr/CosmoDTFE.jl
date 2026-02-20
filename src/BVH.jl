# BVH - Bounding Volume Hierarchy for fast spatial queries
#
# Provides BVHTree abstract type with BVHLeaf and BVHNode, plus
# BoundingVolumeHierarchy wrapper for tetrahedra lookups.


# When operating in BVH, I use simplices rank 3 arrays. This is for memory efficiency and speed.
# Because I realize it can be confusing, I will explain the convention here:
# Simplices have dimensions [3,N,4] where the first dimension is x, y or z,
# the second is the index of the simplex and the third is the vertex of the simplex.
# It may be odd that the vertices are in the third dimension and the second counts up, but that
# way converting to bounding boxes is more simple.

"""
    BVHTree

Abstract type for BVH tree nodes.

    Leafs contain the indices of the tetrahedra for small memory usage.
    Nodes contain children nodes and depth for fast queries.
"""
abstract type BVHTree end

"""
    BVHLeaf

Leaf node of the BVH tree. Contains the indices of the tetrahedra for small memory usage.
"""
struct BVHLeaf <: BVHTree
    data::Vector{Int}
end

"""
    BVHNode

Internal node of the BVH tree. Contains children nodes and depth for fast queries.
"""
struct BVHNode <: BVHTree
    depth::Int
    leftChild::BVHTree
    rightChild::BVHTree
end

"""
    BoundingVolumeHierarchy

BVH tree combined with the bounding box for metadata stored only once.
"""
struct BoundingVolumeHierarchy
    tree::BVHTree
    bbox::Matrix{Float64}
end



function cornerSimplexMatrix(simplex)
    return hcat(minimum(simplex, dims=2), maximum(simplex, dims=2))
end


# Equality operators for ease
function Base.:(==)(a::BVHLeaf, b::BVHLeaf)
    return Set(a.data) == Set(b.data)
end

function Base.:(==)(a::BVHNode, b::BVHNode)
    return a.depth == b.depth &&
           a.leftChild == b.leftChild &&
           a.rightChild == b.rightChild
end

function Base.:(==)(a::BoundingVolumeHierarchy, b::BoundingVolumeHierarchy)
    return a.tree == b.tree && a.bbox == b.bbox
end


"""
    getBoxes(data)

Compute the bounding boxes for all simplices in the dataset.

    The boxes follow the size [3,N,2] where the first dimension is x, y or z,
    the second is the index of the simplex and the third is the minimum or maximum value.
"""
function getBoxes(data::Vector)
    permutedims(stack(cornerSimplexMatrix.(data)), (1, 3, 2))
end

function getBoxes(data::AbstractArray{T,3}) where T # casue there is no general Float type

    minVals = dropdims(minimum(data, dims=3), dims=3)
    maxVals = dropdims(maximum(data, dims=3), dims=3)
    return cat(minVals, maxVals, dims=3)
end


"""
    BoundingVolumeHierarchy(data, depth; box=nothing)

Construct a BVH from simplex data.
"""
function BoundingVolumeHierarchy(data, depth::Int; box=nothing)
    boxes = getBoxes(data)

    if box === nothing
        allMins = @view boxes[:, :, 1]
        allMaxs = @view boxes[:, :, 2]

        globalMin = minimum(allMins, dims=2)
        globalMax = maximum(allMaxs, dims=2)

        box = hcat(globalMin, globalMax)
    end

    tree = generateBvhTree(boxes, depth, box)
    return BoundingVolumeHierarchy(tree, box)
end


"""
    generateBvhTree(boxes, depth::Int, limBox::Matrix) - general case

Generates tree based on number of boxes (total number of tetrahedra).
Boxes are used to simplify logic and speedup.
"""
function generateBvhTree(boxes, depth::Int, limBox::Matrix)
    indices = 1:size(boxes, 2)
    return generateBvhTree(boxes, depth, limBox, indices)
end

"""
    generateBvhTree(boxes, depth::Int, limBox::Matrix, indices) - main case

Generates BVH tree recursively. If depth is 0 or number of tetrahedra is less than 2, returns leaf node.
Otherwise, splits tetrahedra into left and right children and recursively generates left and right children.
"""
function generateBvhTree(boxes, depth::Int, limBox::Matrix, indices)
    if depth == 0 || size(indices, 1) < 2
        return BVHLeaf(collect(indices))
    end

    ax = depth % 3 + 1

    mins = boxes[ax, indices, 1]
    maxs = boxes[ax, indices, 2]

    line = (limBox[ax, 2] + limBox[ax, 1]) / 2

    leftIds = indices[mins.≤line]
    rightIds = indices[maxs.≥line]

    leftBox = copy(limBox)
    leftBox[ax, 2] = line

    rightBox = copy(limBox)
    rightBox[ax, 1] = line

    return BVHNode(depth,
        generateBvhTree(boxes, depth - 1, leftBox, leftIds),
        generateBvhTree(boxes, depth - 1, rightBox, rightIds))
end


# Pretty printing vibe coded
function Base.show(io::IO, ::MIME"text/plain", bvh::BoundingVolumeHierarchy)
    println(io, "BoundingVolumeHierarchy")
    println(io, "  Bounds: ", bvh.bbox)
    print(io, "  Tree:")
    showTree(io, bvh.tree, "  ", true)
end

function showTree(io::IO, node::BVHNode, prefix::String, isLast::Bool)
    println(io)
    print(io, prefix, isLast ? "└─ " : "├─ ")
    print(io, "Node (Depth $(node.depth))")

    newPrefix = prefix * (isLast ? "   " : "│  ")
    showTree(io, node.leftChild, newPrefix, false)
    showTree(io, node.rightChild, newPrefix, true)
end

function showTree(io::IO, leaf::BVHLeaf, prefix::String, isLast::Bool)
    println(io)
    print(io, prefix, isLast ? "└─ " : "├─ ")
    if length(leaf.data) <= 10
        print(io, "Leaf: ", leaf.data)
    else
        print(io, "Leaf ($(length(leaf.data)) elements): ", leaf.data[1:5], "...")
    end
end
