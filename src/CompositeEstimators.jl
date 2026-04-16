# CompositeEstimators.jl
# Implementation of Composite Estimators via BVH for large datasets

abstract type CompositeBVHTree{T} end

"""
    CompositeBVHLeaf

Leaf node of the composite BVH tree containing the actual abstract estimator and its bounding box.
"""
struct CompositeBVHLeaf{T,EstType<:AbstractEstimator{T}} <: CompositeBVHTree{T}
    estimator::EstType
    bounds::Matrix{Float64} # 3x2 matrix [mins maxs]
end

"""
    CompositeBVHNode

Internal node of the composite BVH tree defining the spatial split rule.
Split each axis by the center of mass of the points in the current node.
"""
struct CompositeBVHNode{T} <: CompositeBVHTree{T}
    leftChild::CompositeBVHTree{T}
    rightChild::CompositeBVHTree{T}
    splitAxis::Int
    splitValue::Float64
    bounds::Matrix{Float64}
end

"""
    CompositeEstimator

Composite version of an AbstractEstimator which dynamically routes evaluation logic
to localized estimators using a Bounding Volume Hierarchy to deal with extremely large datasets.
"""
struct CompositeEstimator{T,EstType<:AbstractEstimator{T}} <: AbstractEstimator{T}
    tree::CompositeBVHTree{T}
    maxPoints::Int
    padding::Float64
end

# Functor implementation for single point evaluation
function (est::CompositeEstimator{T,EstType})(point::AbstractVector{<:Real}) where {T,EstType}
    return evaluateComposite(est.tree, point)
end

function evaluateComposite(node::CompositeBVHNode{T}, point::AbstractVector{<:Real}) where T
    if point[node.splitAxis] <= node.splitValue
        return evaluateComposite(node.leftChild, point)
    else
        return evaluateComposite(node.rightChild, point)
    end
end

function evaluateComposite(leaf::CompositeBVHLeaf{T}, point::AbstractVector{<:Real}) where T
    return leaf.estimator(point)
end

"""
    CompositeEstimator(::Type{DensityEstimator}, points, weights; maxPoints=1_000_000, padding=0.1)

Build a multi-level composite DensityEstimator from an extremely large point cloud.
"""
function CompositeEstimator(::Type{DensityEstimator}, points, weights; maxPoints::Int=1_000_000, padding::Float64=0.1, periodic::Bool=false)
    # The output type of DensityEstimator is Float64
    indices = collect(1:length(points))
    if periodic
        boxMins = zeros(Float64, 3)
        boxMaxs = zeros(Float64, 3)
        for d in 1:3
            boxMins[d] = minimum(points[i][d] for i in indices)
            boxMaxs[d] = maximum(points[i][d] for i in indices)
        end
        boxSpans = boxMaxs .- boxMins
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding, boxMins, boxSpans)
    else
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding)
    end
    return CompositeEstimator{Float64,DensityEstimator}(tree, maxPoints, padding)
end

function CompositeEstimator(::Type{DensityEstimator}, points; maxPoints::Int=1_000_000, padding::Float64=0.1, periodic::Bool=false)
    # Default weighting fallback when no weights are passed
    indices = collect(1:length(points))
    weights = ones(Float64, length(points))
    if periodic
        boxMins = zeros(Float64, 3)
        boxMaxs = zeros(Float64, 3)
        for d in 1:3
            boxMins[d] = minimum(points[i][d] for i in indices)
            boxMaxs[d] = maximum(points[i][d] for i in indices)
        end
        boxSpans = boxMaxs .- boxMins
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding, boxMins, boxSpans)
    else
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding)
    end
    return CompositeEstimator{Float64,DensityEstimator}(tree, maxPoints, padding)
end

"""
    CompositeEstimator(::Type{VelocityEstimator}, points, velocities; maxPoints=1_000_000, padding=0.1, periodic=false)

Build a multi-level composite VelocityEstimator from an extremely large point cloud.
"""
function CompositeEstimator(::Type{VelocityEstimator}, points, velocities; maxPoints::Int=1_000_000, padding::Float64=0.1, periodic::Bool=false)
    # The output type of VelocityEstimator is SVector{3, Float64}
    indices = collect(1:length(points))
    if periodic
        boxMins = zeros(Float64, 3)
        boxMaxs = zeros(Float64, 3)
        for d in 1:3
            boxMins[d] = minimum(points[i][d] for i in indices)
            boxMaxs[d] = maximum(points[i][d] for i in indices)
        end
        boxSpans = boxMaxs .- boxMins
        tree = buildCompositeTree(VelocityEstimator, SVector{3,Float64}, points, velocities, indices, maxPoints, padding, boxMins, boxSpans)
    else
        tree = buildCompositeTree(VelocityEstimator, SVector{3,Float64}, points, velocities, indices, maxPoints, padding)
    end
    return CompositeEstimator{SVector{3,Float64},VelocityEstimator}(tree, maxPoints, padding)
end

# Gradient and Velocity aliases for VelocityEstimator
function velocity(est::CompositeEstimator{SVector{3,Float64},VelocityEstimator}, point::AbstractVector{<:Real})
    est(point)
end

function velocityGradient(est::CompositeEstimator{SVector{3,Float64},VelocityEstimator}, point::AbstractVector{<:Real})
    # Determine which leaf estimator covers this point and query gradient there
    return velocityGradientComposite(est.tree, point)
end

function velocityGradientComposite(node::CompositeBVHNode, point::AbstractVector{<:Real})
    if point[node.splitAxis] <= node.splitValue
        return velocityGradientComposite(node.leftChild, point)
    else
        return velocityGradientComposite(node.rightChild, point)
    end
end

function velocityGradientComposite(leaf::CompositeBVHLeaf, point::AbstractVector{<:Real})
    return velocityGradient(leaf.estimator, point)
end

# Recursive builder for BVH. (Non-periodic) Properties are specific for estimator type - weights for DensityEstimator, velocities for VelocityEstimator
function buildCompositeTree(estType::Type, T::Type, globalPoints, globalProperties, indices::Vector{Int}, maxPoints::Int, padding::Float64)
    mins = zeros(Float64, 3)
    maxs = zeros(Float64, 3)

    # Calculate limits for current indices subset
    for d in 1:3
        mins[d] = minimum(globalPoints[i][d] for i in indices)
        maxs[d] = maximum(globalPoints[i][d] for i in indices)
    end
    bounds = hcat(mins, maxs)

    # Base Case: Reach leaf max points
    if length(indices) <= maxPoints
        paddedMins = zeros(3)
        paddedMaxs = zeros(3)
        for d in 1:3
            span = maxs[d] - mins[d]
            paddedMins[d] = mins[d] - padding * span
            paddedMaxs[d] = maxs[d] + padding * span
        end

        paddedIndices = Int[]
        @inbounds for i in 1:length(globalPoints)
            pt = globalPoints[i]
            if pt[1] >= paddedMins[1] && pt[1] <= paddedMaxs[1] &&
               pt[2] >= paddedMins[2] && pt[2] <= paddedMaxs[2] &&
               pt[3] >= paddedMins[3] && pt[3] <= paddedMaxs[3]
                push!(paddedIndices, i)
            end
        end

        leafPoints = globalPoints[paddedIndices]
        if globalProperties !== nothing
            leafProperties = globalProperties[paddedIndices]
            estimator = estType(leafPoints, leafProperties)
        else
            estimator = estType(leafPoints)
        end

        return CompositeBVHLeaf{T,estType}(estimator, bounds)
    end

    # Spatial partitioning by center of mass
    cm = zeros(3)
    for i in indices
        cm .+= globalPoints[i]
    end
    cm ./= length(indices) # uniform mean weights

    spans = maxs .- mins
    splitAxis = argmax(spans)
    splitValue = cm[splitAxis]

    leftIndices = Int[]
    rightIndices = Int[]
    for i in indices
        if globalPoints[i][splitAxis] <= splitValue
            push!(leftIndices, i)
        else
            push!(rightIndices, i)
        end
    end

    # Generate descendants
    leftChild = buildCompositeTree(estType, T, globalPoints, globalProperties, leftIndices, maxPoints, padding)
    rightChild = buildCompositeTree(estType, T, globalPoints, globalProperties, rightIndices, maxPoints, padding)

    return CompositeBVHNode{T}(leftChild, rightChild, splitAxis, splitValue, bounds)
end

# Recursive builder for BVH (Periodic)
function buildCompositeTree(estType::Type, T::Type, globalPoints, globalProperties, indices::Vector{Int}, maxPoints::Int, padding::Float64, boxMins::Vector{Float64}, boxSpans::Vector{Float64})
    mins = zeros(Float64, 3)
    maxs = zeros(Float64, 3)

    # Calculate limits for current indices subset
    for d in 1:3
        mins[d] = minimum(globalPoints[i][d] for i in indices)
        maxs[d] = maximum(globalPoints[i][d] for i in indices)
    end
    bounds = hcat(mins, maxs)

    # Base Case: Reach leaf max points
    if length(indices) <= maxPoints
        paddedMins = zeros(3)
        paddedMaxs = zeros(3)
        for d in 1:3
            span = maxs[d] - mins[d]
            paddedMins[d] = mins[d] - padding * span
            paddedMaxs[d] = maxs[d] + padding * span
        end

        needsWrapCheck = false
        for d in 1:3
            if paddedMins[d] < boxMins[d] || paddedMaxs[d] > (boxMins[d] + boxSpans[d])
                needsWrapCheck = true
                break
            end
        end

        if needsWrapCheck
            leafPoints = typeof(globalPoints[1])[]
            leafProperties = globalProperties !== nothing ? typeof(globalProperties[1])[] : nothing

            center = zeros(3)
            for d in 1:3
                center[d] = (paddedMins[d] + paddedMaxs[d]) / 2.0
            end

            @inbounds for i in 1:length(globalPoints)
                pt = globalPoints[i]

                shiftedX = pt[1]
                distX = shiftedX - center[1]
                shiftedX = center[1] + (distX - boxSpans[1] * round(distX / boxSpans[1]))

                shiftedY = pt[2]
                distY = shiftedY - center[2]
                shiftedY = center[2] + (distY - boxSpans[2] * round(distY / boxSpans[2]))

                shiftedZ = pt[3]
                distZ = shiftedZ - center[3]
                shiftedZ = center[3] + (distZ - boxSpans[3] * round(distZ / boxSpans[3]))

                if shiftedX >= paddedMins[1] && shiftedX <= paddedMaxs[1] &&
                   shiftedY >= paddedMins[2] && shiftedY <= paddedMaxs[2] &&
                   shiftedZ >= paddedMins[3] && shiftedZ <= paddedMaxs[3]

                    push!(leafPoints, SVector{3,Float64}(shiftedX, shiftedY, shiftedZ))
                    if globalProperties !== nothing
                        push!(leafProperties, globalProperties[i])
                    end
                end
            end

            if globalProperties !== nothing
                estimator = estType(leafPoints, leafProperties)
            else
                estimator = estType(leafPoints)
            end
        else
            paddedIndices = Int[]
            @inbounds for i in 1:length(globalPoints)
                pt = globalPoints[i]
                if pt[1] >= paddedMins[1] && pt[1] <= paddedMaxs[1] &&
                   pt[2] >= paddedMins[2] && pt[2] <= paddedMaxs[2] &&
                   pt[3] >= paddedMins[3] && pt[3] <= paddedMaxs[3]
                    push!(paddedIndices, i)
                end
            end

            leafPoints = globalPoints[paddedIndices]
            if globalProperties !== nothing
                leafProperties = globalProperties[paddedIndices]
                estimator = estType(leafPoints, leafProperties)
            else
                estimator = estType(leafPoints)
            end
        end

        return CompositeBVHLeaf{T,estType}(estimator, bounds)
    end

    # Spatial partitioning by center of mass
    cm = zeros(3)
    for i in indices
        cm .+= globalPoints[i]
    end
    cm ./= length(indices) # uniform mean weights

    spans = maxs .- mins
    splitAxis = argmax(spans)
    splitValue = cm[splitAxis]

    leftIndices = Int[]
    rightIndices = Int[]
    for i in indices
        if globalPoints[i][splitAxis] <= splitValue
            push!(leftIndices, i)
        else
            push!(rightIndices, i)
        end
    end

    # Generate descendants
    leftChild = buildCompositeTree(estType, T, globalPoints, globalProperties, leftIndices, maxPoints, padding, boxMins, boxSpans)
    rightChild = buildCompositeTree(estType, T, globalPoints, globalProperties, rightIndices, maxPoints, padding, boxMins, boxSpans)

    return CompositeBVHNode{T}(leftChild, rightChild, splitAxis, splitValue, bounds)
end

# =============================================================================
# Parallel Overloads
# =============================================================================

"""
    CompositeEstimator(::Type{DensityEstimator}, points, weights, nWorkers::Int; maxPoints, padding, periodic)

Parallel constructor. `nWorkers` controls the number of threads used for tree building.
Use `Threads.nthreads()` to utilize all available threads.
"""
function CompositeEstimator(::Type{DensityEstimator}, points, weights, nWorkers::Int; maxPoints::Int=1_000_000, padding::Float64=0.1, periodic::Bool=false)
    indices = collect(1:length(points))
    if periodic
        boxMins = zeros(Float64, 3)
        boxMaxs = zeros(Float64, 3)
        for d in 1:3
            boxMins[d] = minimum(points[i][d] for i in indices)
            boxMaxs[d] = maximum(points[i][d] for i in indices)
        end
        boxSpans = boxMaxs .- boxMins
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding, boxMins, boxSpans, nWorkers)
    else
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding, nWorkers)
    end
    return CompositeEstimator{Float64,DensityEstimator}(tree, maxPoints, padding)
end

function CompositeEstimator(::Type{DensityEstimator}, points, nWorkers::Int; maxPoints::Int=1_000_000, padding::Float64=0.1, periodic::Bool=false)
    indices = collect(1:length(points))
    weights = ones(Float64, length(points))
    if periodic
        boxMins = zeros(Float64, 3)
        boxMaxs = zeros(Float64, 3)
        for d in 1:3
            boxMins[d] = minimum(points[i][d] for i in indices)
            boxMaxs[d] = maximum(points[i][d] for i in indices)
        end
        boxSpans = boxMaxs .- boxMins
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding, boxMins, boxSpans, nWorkers)
    else
        tree = buildCompositeTree(DensityEstimator, Float64, points, weights, indices, maxPoints, padding, nWorkers)
    end
    return CompositeEstimator{Float64,DensityEstimator}(tree, maxPoints, padding)
end

"""
    CompositeEstimator(::Type{VelocityEstimator}, points, velocities, nWorkers::Int; maxPoints, padding, periodic)

Parallel constructor. `nWorkers` controls the number of threads used for tree building.
"""
function CompositeEstimator(::Type{VelocityEstimator}, points, velocities, nWorkers::Int; maxPoints::Int=1_000_000, padding::Float64=0.1, periodic::Bool=false)
    indices = collect(1:length(points))
    if periodic
        boxMins = zeros(Float64, 3)
        boxMaxs = zeros(Float64, 3)
        for d in 1:3
            boxMins[d] = minimum(points[i][d] for i in indices)
            boxMaxs[d] = maximum(points[i][d] for i in indices)
        end
        boxSpans = boxMaxs .- boxMins
        tree = buildCompositeTree(VelocityEstimator, SVector{3,Float64}, points, velocities, indices, maxPoints, padding, boxMins, boxSpans, nWorkers)
    else
        tree = buildCompositeTree(VelocityEstimator, SVector{3,Float64}, points, velocities, indices, maxPoints, padding, nWorkers)
    end
    return CompositeEstimator{SVector{3,Float64},VelocityEstimator}(tree, maxPoints, padding)
end

# Parallel builder for BVH (Non-periodic)
function buildCompositeTree(estType::Type, T::Type, globalPoints, globalProperties, indices::Vector{Int}, maxPoints::Int, padding::Float64, nWorkers::Int)
    # Delegate to serial when budget exhausted or at leaf
    if nWorkers <= 1 || length(indices) <= maxPoints
        return buildCompositeTree(estType, T, globalPoints, globalProperties, indices, maxPoints, padding)
    end

    mins = zeros(Float64, 3)
    maxs = zeros(Float64, 3)
    for d in 1:3
        mins[d] = minimum(globalPoints[i][d] for i in indices)
        maxs[d] = maximum(globalPoints[i][d] for i in indices)
    end
    bounds = hcat(mins, maxs)

    # Spatial partitioning by center of mass
    cm = zeros(3)
    for i in indices
        cm .+= globalPoints[i]
    end
    cm ./= length(indices)

    spans = maxs .- mins
    splitAxis = argmax(spans)
    splitValue = cm[splitAxis]

    leftIndices = Int[]
    rightIndices = Int[]
    for i in indices
        if globalPoints[i][splitAxis] <= splitValue
            push!(leftIndices, i)
        else
            push!(rightIndices, i)
        end
    end

    # Fork-join: spawn right child, compute left on current thread
    nRight = nWorkers ÷ 2
    nLeft = nWorkers - nRight
    rightFuture = Threads.@spawn buildCompositeTree(estType, T, globalPoints, globalProperties, rightIndices, maxPoints, padding, nRight)
    leftChild = buildCompositeTree(estType, T, globalPoints, globalProperties, leftIndices, maxPoints, padding, nLeft)
    rightChild = fetch(rightFuture)

    return CompositeBVHNode{T}(leftChild, rightChild, splitAxis, splitValue, bounds)
end

# Parallel builder for BVH (Periodic)
function buildCompositeTree(estType::Type, T::Type, globalPoints, globalProperties, indices::Vector{Int}, maxPoints::Int, padding::Float64, boxMins::Vector{Float64}, boxSpans::Vector{Float64}, nWorkers::Int)
    # Delegate to serial when budget exhausted or at leaf
    if nWorkers <= 1 || length(indices) <= maxPoints
        return buildCompositeTree(estType, T, globalPoints, globalProperties, indices, maxPoints, padding, boxMins, boxSpans)
    end

    mins = zeros(Float64, 3)
    maxs = zeros(Float64, 3)
    for d in 1:3
        mins[d] = minimum(globalPoints[i][d] for i in indices)
        maxs[d] = maximum(globalPoints[i][d] for i in indices)
    end
    bounds = hcat(mins, maxs)

    # Spatial partitioning by center of mass
    cm = zeros(3)
    for i in indices
        cm .+= globalPoints[i]
    end
    cm ./= length(indices)

    spans = maxs .- mins
    splitAxis = argmax(spans)
    splitValue = cm[splitAxis]

    leftIndices = Int[]
    rightIndices = Int[]
    for i in indices
        if globalPoints[i][splitAxis] <= splitValue
            push!(leftIndices, i)
        else
            push!(rightIndices, i)
        end
    end

    # Fork-join: spawn right child, compute left on current thread
    nRight = nWorkers ÷ 2
    nLeft = nWorkers - nRight
    rightFuture = Threads.@spawn buildCompositeTree(estType, T, globalPoints, globalProperties, rightIndices, maxPoints, padding, boxMins, boxSpans, nRight)
    leftChild = buildCompositeTree(estType, T, globalPoints, globalProperties, leftIndices, maxPoints, padding, boxMins, boxSpans, nLeft)
    rightChild = fetch(rightFuture)

    return CompositeBVHNode{T}(leftChild, rightChild, splitAxis, splitValue, bounds)
end
