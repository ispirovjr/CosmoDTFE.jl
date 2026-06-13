using Test
using StaticArrays
using LinearAlgebra
using CosmoDTFE

@testset "VelocityEstimator" begin
    points = vec([
        SVector{3,Float64}(xValue, yValue, zValue)
        for xValue in 1:5, yValue in 1:5, zValue in 1:5
    ])

    @testset "Constant Velocity Field" begin
        constantVelocity = SVector(1.0, 2.0, 3.0)
        velocities = [constantVelocity for _ in points]

        estimator = VelocityEstimator(points, velocities)
        queryPoint = SVector(2.5, 2.5, 2.5)
        velocityValue, divergence, shear, curl = velocityGradient(estimator, queryPoint)

        @test velocityValue ≈ constantVelocity
        @test divergence ≈ 0.0 atol = 1e-10
        @test norm(shear) ≈ 0.0 atol = 1e-10
        @test norm(curl) ≈ 0.0 atol = 1e-10
    end

    @testset "Pure Expansion" begin
        velocities = [point for point in points]

        estimator = VelocityEstimator(points, velocities)
        queryPoint = SVector(2.5, 2.5, 2.5)
        velocityValue, divergence, shear, curl = velocityGradient(estimator, queryPoint)

        @test velocityValue ≈ queryPoint
        @test divergence ≈ 3.0
        @test norm(shear) ≈ 0.0 atol = 1e-10
        @test norm(curl) ≈ 0.0 atol = 1e-10
    end

    @testset "Solid Body Rotation" begin
        velocities = [SVector(-point[2], point[1], 0.0) for point in points]

        estimator = VelocityEstimator(points, velocities)
        queryPoint = SVector(2.5, 2.5, 2.5)
        velocityValue, divergence, shear, curl = velocityGradient(estimator, queryPoint)

        expectedVelocity = SVector(-queryPoint[2], queryPoint[1], 0.0)
        @test velocityValue ≈ expectedVelocity
        @test divergence ≈ 0.0 atol = 1e-10
        @test curl ≈ SVector(0.0, 0.0, 2.0)
        @test norm(shear) ≈ 0.0 atol = 1e-10
    end

    @testset "Pure Shear" begin
        velocities = [SVector(point[2], point[1], 0.0) for point in points]

        estimator = VelocityEstimator(points, velocities)
        queryPoint = SVector(2.5, 2.5, 2.5)
        velocityValue, divergence, shear, curl = velocityGradient(estimator, queryPoint)

        expectedVelocity = SVector(queryPoint[2], queryPoint[1], 0.0)
        expectedShear = SMatrix{3,3,Float64}(0, 1, 0, 1, 0, 0, 0, 0, 0)

        @test velocityValue ≈ expectedVelocity
        @test divergence ≈ 0.0 atol = 1e-10
        @test norm(curl) ≈ 0.0 atol = 1e-10
        @test shear ≈ expectedShear
    end

    @testset "Bounds Checking" begin
        badPoints = vec([
            SVector{3,Float64}(xValue, yValue, zValue)
            for xValue in 1:2, yValue in 1:2, zValue in 1:2
        ])
        badVelocities = [SVector(0.0, 0.0, 0.0) for _ in 1:7]

        @test_throws ArgumentError VelocityEstimator(badPoints, badVelocities)

        densityEstimator = DensityEstimator(badPoints)
        @test_throws ArgumentError VelocityEstimator(densityEstimator, badVelocities)
    end
end
