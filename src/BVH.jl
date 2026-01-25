# BVH - Bounding Volume Hierarchy for fast spatial queries
#
# Provides BVHTree abstract type with BVHLeaf and BVHNode, plus
# BoundingVolumeHierarchy wrapper for tetrahedra lookups.

abstract type BVHTree end

struct BVHLeaf <: BVHTree
    data::Vector{Int}
end

struct BVHNode <: BVHTree
    depth::Int
    leftChild::BVHTree
    rightChild::BVHTree
end


"""
    BoundingVolumeHierarchy

BVH tree for fast point-in-tetrahedron queries.

# Fields
- `tree::BVHTree`: Root of the BVH tree
- `bbox::Matrix{Float64}`: Bounding box [3×2] matrix
"""
struct BoundingVolumeHierarchy
    tree::BVHTree
    bbox::Matrix{Float64}
end

function generateBvhTree(boxes, depth::Int, limBox::Matrix)
    indices = 1:size(boxes, 2)
    return generateBvhTree(boxes, depth, limBox, indices)
end

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


function BoundingVolumeHierarchy(data::Vector, depth::Int)
    boxes = stack([cornerSimplexMatrix(simplex) for simplex in data])

    minima = (minimum(boxes[1, 1, :]), minimum(boxes[2, 1, :]), minimum(boxes[3, 1, :]))
    maxima = (maximum(boxes[1, 2, :]), maximum(boxes[2, 2, :]), maximum(boxes[3, 2, :]))

    box = stack([minima, maxima])
    tree = generateBvhTree(boxes, depth, box)

    return BoundingVolumeHierarchy(tree, box)
end

function BoundingVolumeHierarchy(data::Vector, depth::Int, box::Matrix)
    boxes = stack([cornerSimplexMatrix(simplex) for simplex in data])

    tree = generateBvhTree(boxes, depth, box)
    return BoundingVolumeHierarchy(tree, box)
end


function BoundingVolumeHierarchy(data::Array, depth::Int, box::Matrix)
    boxes = cat(minimum(data, dims=3), maximum(data, dims=3), dims=3)
    tree = generateBvhTree(boxes, depth, box)
    return BoundingVolumeHierarchy(tree, box)
end

function BoundingVolumeHierarchy(data::Array, depth::Int)
    boxes = cat(minimum(data, dims=3), maximum(data, dims=3), dims=3)
    box = hcat(minimum(boxes[:, :, 1], dims=2), maximum(boxes[:, :, 2], dims=2))
    tree = generateBvhTree(boxes, depth, box)
    return BoundingVolumeHierarchy(tree, box)
end


function cornerSimplexMatrix(simplex)
    return hcat(minimum(simplex, dims=1)', maximum(simplex, dims=1)')
end
