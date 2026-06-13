using Test
using CosmoDTFE

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

    @testset "summed interpolation over multiple scalar streams is exact" begin
        sourcePoints = sourceFixture
        warpedPoints = overlappingWarpedFixture
        values = [10.0, 20.0, 30.0, 40.0, 100.0]
        estimator = PhaseSpaceEstimator(sourcePoints, warpedPoints, values; depth=2)
        queryPoint = Point3(0.25, 0.25, 0.25)

        @test streamNumber(estimator, queryPoint) == 2
        @test estimator(queryPoint) ≈ 65.0
        @test estimator(Point3(2.0, 2.0, 2.0)) == 0.0
    end

    @testset "vector streams are summed generically" begin
        sourcePoints = sourceFixture
        warpedPoints = overlappingWarpedFixture
        values = [
            Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0),
            Point3(1.0, 1.0, 1.0),
            Point3(2.0, 2.0, 2.0),
        ]
        estimator = PhaseSpaceEstimator(sourcePoints, warpedPoints, values; depth=2)

        @test estimator(Point3(0.25, 0.25, 0.25)) ≈ Point3(1.25, 1.25, 1.25)
    end

    @testset "public constructor tessellates source coordinates" begin
        sourcePoints = [
            Point3(0.0, 0.0, 0.0), Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0), Point3(1.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0), Point3(1.0, 0.0, 1.0),
            Point3(0.0, 1.0, 1.0), Point3(1.0, 1.0, 1.0),
        ]
        values = ones(Float64, length(sourcePoints))
        estimator = PhaseSpaceEstimator(sourcePoints, sourcePoints, values; depth=4)
        queryPoint = Point3(0.5, 0.5, 0.5)
        streams = streamNumber(estimator, queryPoint)

        @test streams >= 1
        @test estimator(queryPoint) ≈ streams
    end

    @testset "constructor validates source, warped, and value lengths" begin
        sourcePoints = sourceFixture
        warpedPoints = overlappingWarpedFixture

        @test_throws ArgumentError PhaseSpaceEstimator(sourcePoints[1:4], warpedPoints, ones(length(warpedPoints)))
        @test_throws ArgumentError PhaseSpaceEstimator(sourcePoints, warpedPoints, ones(length(warpedPoints) - 1))
        @test_throws ArgumentError PeriodicPhaseSpaceEstimator(sourcePoints[1:4], warpedPoints, ones(length(warpedPoints)); boxSize=Point3(1.0, 1.0, 1.0))
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
        values = ones(Float64, length(sourcePoints))

        copySourcePts, copyWarpedPts, copyValues = CosmoDTFE.periodicPhaseSpaceCopies(
            sourcePoints,
            warpedPoints,
            values,
            Point3(0.0, 0.0, 0.0),
            Point3(1.0, 1.0, 1.0),
            0.0,
        )
        unwrappedEstimator = PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, values; boxSize=Point3(1.0, 1.0, 1.0), padding=0.0, depth=1)

        @test unwrappedEstimator isa PhaseSpaceEstimator{Float64}
        @test copyValues == values
        @test copySourcePts[1] ≈ Point3(-0.05, 0.10, 0.10)
        @test copyWarpedPts[1] ≈ Point3(0.03, 0.10, 0.10)
        @test streamNumber(unwrappedEstimator, Point3(0.05, 0.12, 0.12)) >= 1
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
        values = ones(Float64, length(sourcePoints))
        framedEstimator = PeriodicPhaseSpaceEstimator(sourcePoints, sourcePoints, values; boxSize=Point3(1.0, 1.0, 1.0), padding=0.10, depth=4)
        noFrameEstimator = PeriodicPhaseSpaceEstimator(sourcePoints, sourcePoints, values; boxSize=Point3(1.0, 1.0, 1.0), padding=0.0, depth=4)
        queryPoint = Point3(0.50, 0.50, 0.50)
        streams = streamNumber(framedEstimator, queryPoint)

        @test length(framedEstimator.triangulation.points) > length(sourcePoints)
        @test any(point -> point[1] > 1.0, framedEstimator.triangulation.points)
        @test streams >= 1
        @test isapprox(framedEstimator(queryPoint), streams; atol=1e-12)
        @test noFrameEstimator(Point3(1.50, 0.50, 0.50)) == 0.0
    end
end
