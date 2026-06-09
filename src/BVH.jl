"""
    BVHTree

Abstract type for BVH tree nodes.
"""
abstract type BVHTree end

"""
    BVHLeaf

Leaf node containing tetrahedron indices.
"""
struct BVHLeaf <: BVHTree
    data::Vector{Int}
end

"""
    BVHNode

Internal BVH node with two children and a remaining recursion depth.
"""
struct BVHNode <: BVHTree
    depth::Int
    leftChild::BVHTree
    rightChild::BVHTree
end

"""
    BoundingVolumeHierarchy

BVH tree plus a global `3 x 2` bounding box stored as `[mins maxs]`.
"""
struct BoundingVolumeHierarchy
    tree::BVHTree
    bbox::Matrix{Float64}
end

"""
    cornerSimplexMatrix(simplex)

Return a `3 x 2` bounding box for one tetrahedron.
"""
function cornerSimplexMatrix(simplex)
    return hcat(minimum(simplex, dims=2), maximum(simplex, dims=2))
end

function Base.:(==)(leftLeaf::BVHLeaf, rightLeaf::BVHLeaf)
    return Set(leftLeaf.data) == Set(rightLeaf.data)
end

function Base.:(==)(leftNode::BVHNode, rightNode::BVHNode)
    return leftNode.depth == rightNode.depth &&
           leftNode.leftChild == rightNode.leftChild &&
           leftNode.rightChild == rightNode.rightChild
end

function Base.:(==)(leftBvh::BoundingVolumeHierarchy, rightBvh::BoundingVolumeHierarchy)
    return leftBvh.tree == rightBvh.tree && leftBvh.bbox == rightBvh.bbox
end

"""
    getBoxes(data)

Compute all simplex bounding boxes as an array with size `3 x nTetrahedra x 2`.
"""
function getBoxes(data::Vector)
    return permutedims(stack(cornerSimplexMatrix.(data)), (1, 3, 2))
end

function getBoxes(data::AbstractArray{T,3}) where T
    minVals = dropdims(minimum(data, dims=3), dims=3)
    maxVals = dropdims(maximum(data, dims=3), dims=3)
    return cat(minVals, maxVals, dims=3)
end

"""
    BoundingVolumeHierarchy(data, depth; box=nothing)

Construct a BVH from simplex coordinate data.
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
    generateBvhTree(boxes, depth, limitBox)

Generate a BVH tree from precomputed simplex boxes.
"""
function generateBvhTree(boxes, depth::Int, limitBox::Matrix)
    indices = collect(1:size(boxes, 2))
    return generateBvhTree(boxes, depth, limitBox, indices)
end

function generateBvhTree(boxes, depth::Int, limBox::Matrix, indices)
    if depth == 0 || length(indices) < 2
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

    return BVHNode(
        depth,
        generateBvhTree(boxes, depth - 1, leftBox, leftIds),
        generateBvhTree(boxes, depth - 1, rightBox, rightIds),
    )
end

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

    childPrefix = prefix * (isLast ? "   " : "|  ")
    showTree(io, node.leftChild, childPrefix, false)
    showTree(io, node.rightChild, childPrefix, true)
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
