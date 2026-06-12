using Test
using StaticArrays
using CosmoDTFE

@testset "PhaseSpaceEstimator" begin
    @testset "Summed interpolation over multiple streams" begin
        sourcePoints = [
            Point3(0.0, 0.0, 0.0),
            Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0),
            Point3(0.0, 0.0, -1.0),
        ]
        warpedPoints = [
            Point3(0.0, 0.0, 0.0),
            Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0),
            Point3(0.0, 0.0, 1.0),
        ]
        values = ones(Float64, length(sourcePoints))
        estimator = PhaseSpaceEstimator(sourcePoints, warpedPoints, values; depth=1)
        queryPoint = Point3(0.25, 0.25, 0.25)

        @test streamNumber(estimator, queryPoint) == 2
        @test estimator(queryPoint) ≈ 2.0
        @test estimator(Point3(2.0, 2.0, 2.0)) == 0.0
    end

    @testset "Vector fields are generic" begin
        sourcePoints = [
            Point3(0.0, 0.0, 0.0),
            Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0),
            Point3(0.0, 0.0, -1.0),
        ]
        warpedPoints = [
            Point3(0.0, 0.0, 0.0),
            Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0),
            Point3(0.0, 0.0, 1.0),
        ]
        values = [Point3(1.0, 2.0, 3.0) for _ in sourcePoints]
        estimator = PhaseSpaceEstimator(sourcePoints, warpedPoints, values; depth=1)

        @test estimator(Point3(0.25, 0.25, 0.25)) ≈ Point3(2.0, 4.0, 6.0)
    end

    @testset "Public constructor tessellates source coordinates" begin
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

    @testset "Periodic wrapper unwraps crossing points" begin
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
        est = PeriodicPhaseSpaceEstimator(sourcePoints, warpedPoints, values; boxSize=Point3(1.0, 1.0, 1.0), padding=0.0, depth=1)

        @test est isa PhaseSpaceEstimator{Float64}
        @test any(p -> p[1] > 1.0, est.triangulation.points)
    end

    @testset "Periodic wrapper frames warped boundaries" begin
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
        est = PeriodicPhaseSpaceEstimator(sourcePoints, sourcePoints, values; boxSize=Point3(1.0, 1.0, 1.0), padding=0.10, depth=4)

        @test length(est.triangulation.points) > length(sourcePoints)
        @test any(p -> p[1] > 1.0, est.triangulation.points)

        queryPoint = Point3(0.50, 0.50, 0.50)
        streams = streamNumber(est, queryPoint)
        @test streams >= 1
        @test isapprox(est(queryPoint), streams; atol=1e-12)

        noFrameEst = PeriodicPhaseSpaceEstimator(sourcePoints, sourcePoints, values; boxSize=Point3(1.0, 1.0, 1.0), padding=0.0, depth=4)
        @test noFrameEst(Point3(1.50, 0.50, 0.50)) == 0.0
    end
end
