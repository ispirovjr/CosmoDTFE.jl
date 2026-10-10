"""
    PhaseSpaceEstimator

Phase-space DTFE (PS-DTFE) density estimator (Feldbrugge 2024). The tessellation is generated from `sourcePoints`, then
every tetrahedron is deformed to `warpedPoints`. Vertex densities are computed
from the deformed star volumes and overlapping warped tetrahedra contribute additively.
"""
struct PhaseSpaceEstimator{T} <: AbstractEstimator{T}
    bvh::BoundingVolumeHierarchy
    triangulation::Triangulation3D{T}
    tetrahedra::Matrix{Int32}
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
    PhaseSpaceEstimator(sourcePoints, warpedPoints, masses; depth=9)
    PhaseSpaceEstimator(sourcePoints, warpedPoints; depth=9)

Build a phase-space DTFE density estimator. The Delaunay connectivity is
computed from `sourcePoints` (Lagrangian / initial positions) and kept fixed
while the vertices are moved to `warpedPoints` (Eulerian / final positions).
The vertex density is evaluated after this deformation.Unit masses are used when `masses` is omitted.
"""
function PhaseSpaceEstimator(sourcePoints, warpedPoints, masses::AbstractVector; depth::Int=9)
    sourcePts = toPoint3Vector(sourcePoints)
    warpedPts = toPoint3Vector(warpedPoints)
    length(sourcePts) == length(warpedPts) || throw(ArgumentError("sourcePoints and warpedPoints must have the same length"))
    length(warpedPts) == length(masses) || throw(ArgumentError("masses must have the same length as warpedPoints"))

    _, tetrahedra = tessellate(sourcePts)
    triangulation = Triangulation3D(warpedPts, tetrahedra, masses)
    bvh = BoundingVolumeHierarchy(warpedPts, tetrahedra, depth)

    return PhaseSpaceEstimator{Float64}(bvh, triangulation, tetrahedra)
end

function PhaseSpaceEstimator(sourcePoints, warpedPoints; depth::Int=9)
    return PhaseSpaceEstimator(sourcePoints, warpedPoints, ones(Float64, length(sourcePoints)); depth=depth)
end

"""
    PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, masses; boxSize, boxMin=Point3(0,0,0), padding=0.05, depth=9)
    PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints; boxSize, boxMin=Point3(0,0,0), padding=0.05, depth=9)

Build a periodic phase-space DTFE density estimator. Boundary copies are
selected from minimum-image-unwrapped `warpedPoints`; each copied warped point
carries the matching source point and mass. The result is a normal
`PhaseSpaceEstimator` whose densities are computed on the framed point set.
Vertices at the outer edge of the padding frame have incomplete stars, so
`padding` should exceed the size of the largest warped star near the box faces.
"""
function PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, masses::AbstractVector; boxSize, boxMin=Point3(0.0, 0.0, 0.0), padding::Real=0.05, depth::Int=9)
    sourcePts = toPoint3Vector(sourcePoints)
    warpedPts = toPoint3Vector(warpedPoints)
    pointMasses = collect(Float64, masses)
    length(sourcePts) == length(warpedPts) || throw(ArgumentError("sourcePoints and warpedPoints must have the same length"))
    length(sourcePts) == length(pointMasses) || throw(ArgumentError("masses must have the same length as sourcePoints"))

    boxMin = Point3(boxMin)
    boxSize = Point3(boxSize)
    copySourcePts, copyWarpedPts, copyMasses = periodicPhaseSpaceCopies(sourcePts, warpedPts, pointMasses, boxMin, boxSize, padding)
    return PhaseSpaceEstimator(copySourcePts, copyWarpedPts, copyMasses; depth=depth)
end

function PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints; kwargs...)
    return PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, ones(Float64, length(sourcePoints)); kwargs...)
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
