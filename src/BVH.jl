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
    bbox::SearchBox3D
end

"""
    cornerSimplexMatrix(simplex)

Return a `3 x 2` bounding box for one `3 x 4` tetrahedron coordinate matrix.
"""
function cornerSimplexMatrix(simplex::AbstractMatrix)
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

Compute all simplex bounding boxes as an array with size
`3 x 2 x nTetrahedra`. 
"""
function getBoxes(data::Vector)
    return stack(cornerSimplexMatrix.(data))
end

function getBoxes(data::AbstractArray{T,3}) where T
    minVals = minimum(data, dims=2)
    maxVals = maximum(data, dims=2)
    return cat(minVals, maxVals, dims=2)
end

"""
    getBoxes(points, tetrahedra)

Compute tetrahedron bounding boxes directly from a point vector and a
`4 x nTetrahedra` connectivity matrix.
"""
function getBoxes(points::AbstractVector, tetrahedra::AbstractMatrix{<:Integer})
    boxes = Array{Float64}(undef, 3, 2, size(tetrahedra, 2))

    @inbounds for tetrahedronId in axes(tetrahedra, 2)
        pointOne = points[tetrahedra[1, tetrahedronId]]
        pointTwo = points[tetrahedra[2, tetrahedronId]]
        pointThree = points[tetrahedra[3, tetrahedronId]]
        pointFour = points[tetrahedra[4, tetrahedronId]]

        for axis in 1:3
            boxes[axis, 1, tetrahedronId] = min(
                pointOne[axis],
                pointTwo[axis],
                pointThree[axis],
                pointFour[axis],
            )
            boxes[axis, 2, tetrahedronId] = max(
                pointOne[axis],
                pointTwo[axis],
                pointThree[axis],
                pointFour[axis],
            )
        end
    end

    return boxes
end

function buildBvh(boxes, depth::Int, box)
    if box === nothing
        allMins = @view boxes[:, 1, :]
        allMaxs = @view boxes[:, 2, :]
        globalMin = minimum(allMins, dims=2)
        globalMax = maximum(allMaxs, dims=2)
        box = hcat(globalMin, globalMax)
    end
    box = SearchBox3D(box)

    tree = generateBvhTree(boxes, depth, box)
    return BoundingVolumeHierarchy(tree, box)
end

"""
    BoundingVolumeHierarchy(data, depth; box=nothing)

Construct a BVH from simplex coordinate data.
"""
function BoundingVolumeHierarchy(data, depth::Int; box=nothing)
    boxes = getBoxes(data)
    return buildBvh(boxes, depth, box)
end

"""
    BoundingVolumeHierarchy(points, tetrahedra, depth; box=nothing)

Construct a BVH directly from point coordinates and a `4 x nTetrahedra`
connectivity matrix.
"""
function BoundingVolumeHierarchy(points::AbstractVector, tetrahedra::AbstractMatrix{<:Integer}, depth::Int; box=nothing)
    boxes = getBoxes(points, tetrahedra)
    return buildBvh(boxes, depth, box)
end

"""
    generateBvhTree(boxes, depth, limitBox)

Generate a BVH tree from precomputed simplex boxes.
"""
function generateBvhTree(boxes, depth::Int, limitBox::SearchBox3D)
    indices = collect(1:size(boxes, 3))
    return generateBvhTree(boxes, depth, limitBox, indices)
end

function generateBvhTree(boxes, depth::Int, limBox::SearchBox3D, indices)
    if depth == 0 || length(indices) < 2
        return BVHLeaf(collect(indices))
    end

    ax = depth % 3 + 1
    mins = boxes[ax, 1, indices]
    maxs = boxes[ax, 2, indices]
    line = (limBox[ax, 2] + limBox[ax, 1]) / 2
    
    leftIds = indices[mins.≤line]
    rightIds = indices[maxs.≥line]

    leftBox = setindex(limBox, line, ax, 2)
    rightBox = setindex(limBox, line, ax, 1)

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
