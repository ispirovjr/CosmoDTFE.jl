using Test
using StaticArrays
using LinearAlgebra
using CosmoDTFE

function makePeriodicGrid(nSide=4)
    axisValues = range(0.05, 0.95, length=nSide)
    return [Point3(xValue, yValue, zValue) for xValue in axisValues for yValue in axisValues for zValue in axisValues]
end

@testset "PeriodicEstimator" begin
    points = makePeriodicGrid()
    weights = ones(Float64, length(points))

    @testset "Density copies boundary points" begin
        est = PeriodicEstimator(DensityEstimator, points, weights; boxSize=Point3(1.0, 1.0, 1.0), padding=0.12, depth=5)

        @test est isa PeriodicEstimator{Float64,DensityEstimator}
        @test length(est.estimator.triangulation.points) > length(points)

        lowerDensity = est(Point3(0.02, 0.50, 0.50))
        wrappedDensity = est(Point3(1.02, 0.50, 0.50))

        @test lowerDensity > 0.0
        @test isfinite(lowerDensity)
        @test lowerDensity ≈ wrappedDensity
    end

    @testset "Unweighted density constructor" begin
        est = PeriodicEstimator(DensityEstimator, points; boxSize=Point3(1.0, 1.0, 1.0), padding=0.12, depth=5)
        densityValue = est(Point3(-0.02, 0.50, 0.50))

        @test densityValue > 0.0
        @test isfinite(densityValue)
    end

    @testset "Velocity wrapping preserves constant field" begin
        constantVelocity = Point3(2.0, -1.0, 0.5)
        velocities = [constantVelocity for _ in points]
        est = PeriodicEstimator(VelocityEstimator, points, velocities; boxSize=Point3(1.0, 1.0, 1.0), padding=0.12, depth=5)

        velocityValue = est(Point3(1.02, 0.50, 0.50))
        velocityInterp, divergence, shear, vorticity = velocityGradient(est, Point3(-0.02, 0.50, 0.50))

        @test velocityValue ≈ constantVelocity
        @test velocityInterp ≈ constantVelocity
        @test divergence ≈ 0.0 atol = 1e-10
        @test norm(shear) ≈ 0.0 atol = 1e-10
        @test norm(vorticity) ≈ 0.0 atol = 1e-10
    end
end
