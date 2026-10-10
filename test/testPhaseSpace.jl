using Test
using CosmoDTFE

# Integral of the piecewise-linear DTFE field over all tetrahedra:
# ∫ ρ(x) dx = ∑_T V(T) * (1/4) * ∑_{i ∈ T} rhoStar[i]
function integratedMass(est::PhaseSpaceEstimator)
    total = 0.0
    for tet in eachcol(est.tetrahedra)
        volume = CosmoDTFE.computeVolume(est.triangulation.points[tet])
        total += volume * sum(est.triangulation.rhoStar[tet]) / 4
    end
    return total
end

@testset "PhaseSpaceEstimator" begin
    sourceFixture = [
        Point3(0.0, 0.0, 0.0),
        Point3(1.0, 0.0, 0.0),
        Point3(0.0, 1.0, 0.0),
        Point3(0.0, 0.0, 1.0),
        Point3(0.0, 0.0, -1.0),
    ]

    overlappingWarpedFixture = [
        Point3(0.0, 0.0, 0.0),
        Point3(1.0, 0.0, 0.0),
        Point3(0.0, 1.0, 0.0),
        Point3(0.0, 0.0, 1.0),
        Point3(0.0, 0.0, 1.0),
    ]

    jitteredGrid = vec([
        Point3(i + 1e-3 * sin(7i + 3j + k), j + 1e-3 * cos(5i + j + 2k), k + 1e-3 * sin(i + 11j + 13k))
        for i in 0:4, j in 0:4, k in 0:4
    ])
    gridMasses = [1.0 + 0.1 * mod(pointId, 3) for pointId in eachindex(jitteredGrid)]

    @testset "estimator calculates physical densities rather than raw value interpolation" begin
        # Vertices 1, 2, 3 belong to 2 tetrahedra (star volume = 2 * 1/6 = 1/3)
        # Vertices 4, 5 belong to 1 tetrahedron (star volume = 1/6)
        masses = [10.0, 20.0, 30.0, 40.0, 100.0]
        estimator = PhaseSpaceEstimator(sourceFixture, overlappingWarpedFixture, masses; depth=2)
        queryPoint = Point3(0.25, 0.25, 0.25)

        # 1) rhoStar must reflect 4 * mass / V(W_i), NOT the raw values
        expectedRhoStar = [
            4 * 10.0 / (1/3),   # 120.0
            4 * 20.0 / (1/3),   # 240.0
            4 * 30.0 / (1/3),   # 360.0
            4 * 40.0 / (1/6),   # 960.0
            4 * 100.0 / (1/6),  # 2400.0
        ]
        @test estimator.triangulation.rhoStar ≈ expectedRhoStar
        @test estimator.triangulation.rhoStar != masses

        # 2) Density at queryPoint must be the sum over both streams:
        # Stream 1: 0.25 * (120 + 240 + 360 + 960) = 420.0
        # Stream 2: 0.25 * (120 + 240 + 360 + 2400) = 780.0
        # Total density = 1200.0.
        # (In contrast, raw value interpolation would unphysically yield 65.0)
        @test streamNumber(estimator, queryPoint) == 2
        @test estimator(queryPoint) ≈ 1200.0
        @test estimator(queryPoint) != 65.0
        @test estimator(Point3(2.0, 2.0, 2.0)) == 0.0

        # 3) Total mass conservation holds because rhoStar uses volume of the umbrella/star
        @test integratedMass(estimator) ≈ sum(masses)
    end

    @testset "density scales inversely with volume under deformation" begin
        scale = 2.0
        warpedGrid = [scale * point for point in jitteredGrid]

        undeformedEstimator = PhaseSpaceEstimator(jitteredGrid, jitteredGrid, gridMasses; depth=4)
        expandedEstimator = PhaseSpaceEstimator(jitteredGrid, warpedGrid, gridMasses; depth=4)

        # Under homogeneous 3D expansion by factor s, volumes scale as s^3 = 8
        # Therefore, density must drop by exactly factor s^3 = 8
        @test expandedEstimator.triangulation.rhoStar ≈ undeformedEstimator.triangulation.rhoStar ./ (scale^3)

        queryPoint = Point3(1.5, 1.7, 1.6)
        @test expandedEstimator(scale * queryPoint) ≈ undeformedEstimator(queryPoint) / (scale^3)
        # If it were raw value interpolation, the value at (scale * queryPoint) would equal queryPoint:
        @test !(expandedEstimator(scale * queryPoint) ≈ undeformedEstimator(queryPoint))
    end

    @testset "default constructor uses unit masses" begin
        unitEstimator = PhaseSpaceEstimator(sourceFixture, overlappingWarpedFixture; depth=2)
        explicitEstimator = PhaseSpaceEstimator(sourceFixture, overlappingWarpedFixture, ones(Float64, length(sourceFixture)); depth=2)

        @test unitEstimator.triangulation.rhoStar ≈ explicitEstimator.triangulation.rhoStar
        @test integratedMass(unitEstimator) ≈ length(sourceFixture)
    end

    @testset "masses scale vertex densities linearly" begin
        masses = [1.0, 2.0, 3.0, 4.0, 5.0]
        unitEstimator = PhaseSpaceEstimator(sourceFixture, overlappingWarpedFixture; depth=2)
        massEstimator = PhaseSpaceEstimator(sourceFixture, overlappingWarpedFixture, masses; depth=2)

        @test massEstimator.triangulation.rhoStar ≈ masses .* unitEstimator.triangulation.rhoStar
    end

    @testset "constructor validates input lengths" begin
        @test_throws ArgumentError PhaseSpaceEstimator(sourceFixture[1:4], overlappingWarpedFixture)
        @test_throws ArgumentError PhaseSpaceEstimator(sourceFixture, overlappingWarpedFixture, ones(4))
        @test_throws ArgumentError PeriodicPhaseSpaceEstimator(sourceFixture[1:4], overlappingWarpedFixture; boxSize=Point3(1.0, 1.0, 1.0))
        @test_throws ArgumentError PeriodicPhaseSpaceEstimator(sourceFixture, overlappingWarpedFixture, ones(4); boxSize=Point3(1.0, 1.0, 1.0))
    end

    @testset "periodic phase-space unwraps and reframes crossing warped points" begin
        sourcePoints = [
            Point3(0.95, 0.10, 0.10),
            Point3(0.20, 0.85, 0.10),
            Point3(0.20, 0.10, 0.85),
            Point3(0.20, 0.10, 0.10),
        ]
        warpedPoints = [
            Point3(0.03, 0.10, 0.10),
            sourcePoints[2],
            sourcePoints[3],
            sourcePoints[4],
        ]
        masses = ones(Float64, length(sourcePoints))

        unwrappedEstimator = PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, masses; boxSize=Point3(1.0, 1.0, 1.0), padding=0.0, depth=1)

        @test unwrappedEstimator isa PhaseSpaceEstimator{Float64}
        @test streamNumber(unwrappedEstimator, Point3(0.05, 0.12, 0.12)) >= 1

        # Single tetrahedron: rhoStar = 4 * m / V
        copySourcePts, copyWarpedPts, _ = CosmoDTFE.periodicPhaseSpaceCopies(
            sourcePoints, warpedPoints, masses, Point3(0.0, 0.0, 0.0), Point3(1.0, 1.0, 1.0), 0.0
        )
        tetVolume = CosmoDTFE.computeVolume(copyWarpedPts)
        @test unwrappedEstimator.triangulation.rhoStar ≈ fill(4.0 / tetVolume, 4)
    end

    @testset "periodic phase-space frames warped boundaries" begin
        sourcePoints = [
            Point3(0.04, 0.20, 0.20),
            Point3(0.80, 0.20, 0.20),
            Point3(0.20, 0.80, 0.20),
            Point3(0.20, 0.20, 0.80),
            Point3(0.80, 0.80, 0.20),
            Point3(0.80, 0.20, 0.80),
            Point3(0.20, 0.80, 0.80),
            Point3(0.80, 0.80, 0.80),
        ]
        masses = collect(1.0:length(sourcePoints))
        boxSize = Point3(1.0, 1.0, 1.0)
        framedEstimator = PeriodicPhaseSpaceEstimator(sourcePoints, sourcePoints, masses; boxSize=boxSize, padding=0.10, depth=4)
        noFrameEstimator = PeriodicPhaseSpaceEstimator(sourcePoints, sourcePoints, masses; boxSize=boxSize, padding=0.0, depth=4)
        queryPoint = Point3(0.52, 0.52, 0.52)

        copySourcePts, copyWarpedPts, copyMasses = CosmoDTFE.periodicPhaseSpaceCopies(
            sourcePoints, sourcePoints, masses, Point3(0.0, 0.0, 0.0), boxSize, 0.10,
        )
        referenceEstimator = PhaseSpaceEstimator(copySourcePts, copyWarpedPts, copyMasses; depth=4)

        @test length(framedEstimator.triangulation.points) > length(sourcePoints)
        @test any(point -> point[1] > 1.0, framedEstimator.triangulation.points)
        @test streamNumber(framedEstimator, queryPoint) >= 1
        @test framedEstimator.triangulation.rhoStar ≈ referenceEstimator.triangulation.rhoStar
        @test framedEstimator(queryPoint) ≈ referenceEstimator(queryPoint)
        @test framedEstimator(queryPoint) > 0.0
        @test noFrameEstimator(Point3(1.50, 0.50, 0.50)) == 0.0
    end
end
