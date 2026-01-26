# BVH - Bounding Volume Hierarchy for fast spatial queries
#
# Provides BVHTree abstract type with BVHLeaf and BVHNode, plus
# BoundingVolumeHierarchy wrapper for tetrahedra lookups.

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
    return hcat(minimum(simplex, dims=1)', maximum(simplex, dims=1)')
end


"""
    getBoxes(data)

Compute the bounding boxes for all simplices in the dataset.
"""
function getBoxes(data::Vector)
    permutedims(stack(cornerSimplexMatrix.(data)), (1, 3, 2))
end

function getBoxes(data::AbstractArray{T,3}) where T
    cat(minimum(data, dims=3), maximum(data, dims=3), dims=3)
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
