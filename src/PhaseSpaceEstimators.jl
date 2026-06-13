"""
    PhaseSpaceEstimator

Generic DTFE estimator for warped tessellations. The tessellation is generated
from `sourcePoints`, then each tetrahedron is evaluated in `warpedPoints`.
Overlapping warped tetrahedra contribute additively.
"""
struct PhaseSpaceEstimator{T} <: AbstractEstimator{T}
    bvh::BoundingVolumeHierarchy
    triangulation::Triangulation3D{T}
    tetrahedra::Matrix{Int}
end

function unwrapWarpedPoint(q::Point3, x::Point3, boxSize::Point3)
    return Point3(
        q[1] + mod(x[1] - q[1] + boxSize[1] / 2, boxSize[1]) - boxSize[1] / 2,
        q[2] + mod(x[2] - q[2] + boxSize[2] / 2, boxSize[2]) - boxSize[2] / 2,
        q[3] + mod(x[3] - q[3] + boxSize[3] / 2, boxSize[3]) - boxSize[3] / 2,
    )
end

function periodicPhaseSpaceCopies(sourcePoints::Vector{Point3}, warpedPoints::Vector{Point3}, values, boxMin::Point3, boxSize::Point3, padding::Real)
    pad = padding * boxSize
    boxMax = boxMin + boxSize
    copySourcePts = Point3[]
    copyWarpedPts = Point3[]
    copyVals = Vector{eltype(values)}()
    sizehint!(copySourcePts, length(sourcePoints))
    sizehint!(copyWarpedPts, length(warpedPoints))
    sizehint!(copyVals, length(values))

    @inbounds for pointId in eachindex(sourcePoints)
        q = sourcePoints[pointId]
        x = unwrapWarpedPoint(q, warpedPoints[pointId], boxSize)

        shift = Point3(0.0, 0.0, 0.0) #in case unwrapping moves x out of box
        x[1] < boxMin[1] && (shift += Point3(boxSize[1], 0.0, 0.0))
        x[1] > boxMax[1] && (shift -= Point3(boxSize[1], 0.0, 0.0))
        x[2] < boxMin[2] && (shift += Point3(0.0, boxSize[2], 0.0))
        x[2] > boxMax[2] && (shift -= Point3(0.0, boxSize[2], 0.0))
        x[3] < boxMin[3] && (shift += Point3(0.0, 0.0, boxSize[3]))
        x[3] > boxMax[3] && (shift -= Point3(0.0, 0.0, boxSize[3]))
        q += shift
        x += shift

        xShifts = [0.0]
        yShifts = [0.0]
        zShifts = [0.0]

        x[1] < boxMin[1] + pad[1] && push!(xShifts, boxSize[1])
        x[1] > boxMax[1] - pad[1] && push!(xShifts, -boxSize[1])

        x[2] < boxMin[2] + pad[2] && push!(yShifts, boxSize[2])
        x[2] > boxMax[2] - pad[2] && push!(yShifts, -boxSize[2])

        x[3] < boxMin[3] + pad[3] && push!(zShifts, boxSize[3])
        x[3] > boxMax[3] - pad[3] && push!(zShifts, -boxSize[3])

        for xShift in xShifts, yShift in yShifts, zShift in zShifts
            shift = Point3(xShift, yShift, zShift)
            push!(copySourcePts, q + shift)
            push!(copyWarpedPts, x + shift)
            push!(copyVals, values[pointId])
        end
    end

    return copySourcePts, copyWarpedPts, copyVals
end

"""
    PhaseSpaceEstimator(sourcePoints, warpedPoints, values; depth=9)

Build a generic phase-space estimator. The source points define the
tessellation, while warped points define where each tetrahedron is searched and
interpolated.
"""
function PhaseSpaceEstimator(sourcePoints, warpedPoints, values; depth::Int=9)
    sourcePts = toPoint3Vector(sourcePoints)
    warpedPts = toPoint3Vector(warpedPoints)
    interpValues = collect(values)
    length(sourcePts) == length(warpedPts) || throw(ArgumentError("sourcePoints and warpedPoints must have the same length"))
    length(warpedPts) == length(interpValues) || throw(ArgumentError("values must have the same length as warpedPoints"))

    _, tetrahedra = tessellate(sourcePts)
    triangulation = Triangulation3D{eltype(interpValues)}(warpedPts, interpValues)
    warpedCoords = pointsToColumns(warpedPts)
    bvh = BoundingVolumeHierarchy(warpedCoords[:, tetrahedra], depth)

    return PhaseSpaceEstimator{eltype(triangulation.rhoStar)}(
        bvh,
        triangulation,
        Matrix{Int}(tetrahedra),
    )
end

"""
    PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, values; boxSize, boxMin=Point3(0,0,0), padding=0.05, depth=9)

Build a periodic phase-space estimator. Boundary copies are selected from
minimum-image-unwrapped `warpedPoints`; each copied warped point carries the
matching source point and value. The result is a normal `PhaseSpaceEstimator`.
"""
function PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, values; boxSize, boxMin=Point3(0.0, 0.0, 0.0), padding::Real=0.05, depth::Int=9)
    sourcePts = toPoint3Vector(sourcePoints)
    warpedPts = toPoint3Vector(warpedPoints)
    interpValues = collect(values)
    length(sourcePts) == length(warpedPts) || throw(ArgumentError("sourcePoints and warpedPoints must have the same length"))
    length(sourcePts) == length(interpValues) || throw(ArgumentError("values must have the same length as sourcePoints"))

    boxMin = Point3(boxMin)
    boxSize = Point3(boxSize)
    copySourcePts, copyWarpedPts, copyValues = periodicPhaseSpaceCopies(sourcePts, warpedPts, interpValues, boxMin, boxSize, padding)
    return PhaseSpaceEstimator(copySourcePts, copyWarpedPts, copyValues; depth=depth)
end

function (est::PhaseSpaceEstimator{T})(point::AbstractVector{<:Real}) where T
    tetIds = findAllIds(point, est.triangulation.points, est.tetrahedra, est.bvh)
    isempty(tetIds) && return zero(T)

    result = zero(T)
    for tetId in tetIds
        verts = simplexData(est.triangulation.points, est.tetrahedra, tetId)
        interpValues = simplexData(est.triangulation.rhoStar, est.tetrahedra, tetId)
        lambda = barycentricWeights(point, verts...)
        result += interpolateSimplex(lambda, interpValues...)
    end

    return result
end

"""
    streamNumber(estimator, point)

Return the number of warped tetrahedra containing `point`.
"""
function streamNumber(est::PhaseSpaceEstimator, point::AbstractVector{<:Real})
    tetIds = findAllIds(point, est.triangulation.points, est.tetrahedra, est.bvh)
    return length(tetIds)
end
