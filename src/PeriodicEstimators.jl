"""
    PeriodicEstimator

Wrapper for density and velocity estimators on a periodic rectangular box.
The wrapped estimator is built from the original points plus translated copies
near periodic faces.
"""
struct PeriodicEstimator{T,EstimatorType<:AbstractEstimator{T}} <: AbstractEstimator{T}
    estimator::EstimatorType
    boxMin::Point3
    boxSize::Point3
    padding::Float64
end

@inline function wrapPeriodicPoint(point::AbstractVector{<:Real}, boxMin::Point3, boxSize::Point3)
    p = Point3(point)
    return Point3(
        boxMin[1] + mod(p[1] - boxMin[1], boxSize[1]),
        boxMin[2] + mod(p[2] - boxMin[2], boxSize[2]),
        boxMin[3] + mod(p[3] - boxMin[3], boxSize[3]),
    )
end

function periodicCopies(points, props, boxMin::Point3, boxSize::Point3, padding::Real)
    pts = toPoint3Vector(points)
    length(pts) == length(props) || throw(ArgumentError("length of properties must match length of points"))

    pad = padding * boxSize
    copyPts = Point3[]
    copyProps = Vector{eltype(props)}()
    sizehint!(copyPts, length(pts))
    sizehint!(copyProps, length(props))

    @inbounds for pointId in eachindex(pts)
        p = pts[pointId]

        xShifts = [0.0]
        yShifts = [0.0]
        zShifts = [0.0]

        p[1] - boxMin[1] <= pad[1] && push!(xShifts, boxSize[1])
        boxMin[1] + boxSize[1] - p[1] <= pad[1] && push!(xShifts, -boxSize[1])

        p[2] - boxMin[2] <= pad[2] && push!(yShifts, boxSize[2])
        boxMin[2] + boxSize[2] - p[2] <= pad[2] && push!(yShifts, -boxSize[2])

        p[3] - boxMin[3] <= pad[3] && push!(zShifts, boxSize[3])
        boxMin[3] + boxSize[3] - p[3] <= pad[3] && push!(zShifts, -boxSize[3])

        for xShift in xShifts, yShift in yShifts, zShift in zShifts
            push!(copyPts, p + Point3(xShift, yShift, zShift))
            push!(copyProps, props[pointId])
        end
    end

    return copyPts, copyProps
end

"""
    PeriodicEstimator(DensityEstimator, points, weights; boxSize, boxMin=Point3(0,0,0), padding=0.0, depth=9)
    PeriodicEstimator(DensityEstimator, points; boxSize, boxMin=Point3(0,0,0), padding=0.0, depth=9)

Build a periodic density estimator.
"""
function PeriodicEstimator(::Type{DensityEstimator}, points, weights::AbstractVector; boxSize, boxMin=Point3(0.0, 0.0, 0.0), padding::Real=0.0, depth::Int=9)
    boxMin = Point3(boxMin)
    boxSize = Point3(boxSize)
    copyPts, copyWeights = periodicCopies(points, weights, boxMin, boxSize, padding)
    est = DensityEstimator(copyPts, copyWeights; depth=depth)
    return PeriodicEstimator{Float64,DensityEstimator}(est, boxMin, boxSize, Float64(padding))
end

function PeriodicEstimator(::Type{DensityEstimator}, points; boxSize, boxMin=Point3(0.0, 0.0, 0.0), padding::Real=0.0, depth::Int=9)
    weights = ones(Float64, length(points))
    return PeriodicEstimator(DensityEstimator, points, weights; boxSize=boxSize, boxMin=boxMin, padding=padding, depth=depth)
end

"""
    PeriodicEstimator(VelocityEstimator, points, velocities; boxSize, boxMin=Point3(0,0,0), padding=0.0, depth=9)

Build a periodic velocity estimator.
"""
function PeriodicEstimator(::Type{VelocityEstimator}, points, velocities::AbstractVector; boxSize, boxMin=Point3(0.0, 0.0, 0.0), padding::Real=0.0, depth::Int=9)
    boxMin = Point3(boxMin)
    boxSize = Point3(boxSize)
    copyPts, copyVels = periodicCopies(points, velocities, boxMin, boxSize, padding)
    est = VelocityEstimator(copyPts, copyVels; depth=depth)
    return PeriodicEstimator{SVector{3,Float64},VelocityEstimator}(est, boxMin, boxSize, Float64(padding))
end

function (est::PeriodicEstimator{T})(point::AbstractVector{<:Real}) where T
    wrappedPoint = wrapPeriodicPoint(point, est.boxMin, est.boxSize)
    return est.estimator(wrappedPoint)
end

velocity(est::PeriodicEstimator{SVector{3,Float64},VelocityEstimator}, point::AbstractVector{<:Real}) = est(point)

function velocityGradient(est::PeriodicEstimator{SVector{3,Float64},VelocityEstimator}, point::AbstractVector{<:Real})
    wrappedPoint = wrapPeriodicPoint(point, est.boxMin, est.boxSize)
    return velocityGradient(est.estimator, wrappedPoint)
end
