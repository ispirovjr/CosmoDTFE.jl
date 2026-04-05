using Test
using StaticArrays
using LinearAlgebra
using JuliaDTFE

# Setup a simple 3D grid of points to tessellate
function make_grid_points(N=5)
    points = [SVector{3,Float64}(x, y, z) for x in 1:N, y in 1:N, z in 1:N]
    return vec(points)
end

@testset "VelocityEstimator Tests" begin
    points = make_grid_points(5) # 5x5x5 grid

    @testset "Constant Velocity Field" begin
        # v(x) = (1, 2, 3)
        const_vel = SVector(1.0, 2.0, 3.0)
        velocities = [const_vel for _ in points]

        est = VelocityEstimator(points, velocities)

        # Test inside a known bounds (center of grid)
        query_pt = SVector(2.5, 2.5, 2.5)

        v, div, shear, curl = velocityGradient(est, query_pt)

        @test v ≈ const_vel
        @test div ≈ 0.0 atol = 1e-10
        @test norm(shear) ≈ 0.0 atol = 1e-10
        @test norm(curl) ≈ 0.0 atol = 1e-10
    end

    @testset "Pure Expansion (Linear Field)" begin
        # v(x) = x
        velocities = [p for p in points]

        est = VelocityEstimator(points, velocities)
        query_pt = SVector(2.5, 2.5, 2.5)

        v, div, shear, curl = velocityGradient(est, query_pt)

        # v(x) should be x
        @test v ≈ query_pt

        # Div v = 3
        @test div ≈ 3.0

        # Shear should be 0 because J = I, so S = I, and σ = I - (3/3)I = 0
        @test norm(shear) ≈ 0.0 atol = 1e-10

        # Curl should be 0
        @test norm(curl) ≈ 0.0 atol = 1e-10
    end

    @testset "Solid Body Rotation" begin
        # ω = (0, 0, 1)
        # v = ω × x = (-y, x, 0)
        velocities = [SVector(-p[2], p[1], 0.0) for p in points]

        est = VelocityEstimator(points, velocities)
        query_pt = SVector(2.5, 2.5, 2.5) # Center

        v, div, shear, curl = velocityGradient(est, query_pt)

        expected_v = SVector(-query_pt[2], query_pt[1], 0.0)
        @test v ≈ expected_v

        # Div = 0
        @test div ≈ 0.0 atol = 1e-10

        # Curl = 2ω = (0, 0, 2)
        @test curl ≈ SVector(0.0, 0.0, 2.0)

        # J = [0 -1 0; 1 0 0; 0 0 0]
        # S = 0.5(J + J') = 0
        # Shear = 0
        @test norm(shear) ≈ 0.0 atol = 1e-10
    end

    @testset "Pure Shear" begin
        # v(x,y,z) = (y, x, 0)
        velocities = [SVector(p[2], p[1], 0.0) for p in points]

        est = VelocityEstimator(points, velocities)
        query_pt = SVector(2.5, 2.5, 2.5)

        v, div, shear, curl = velocityGradient(est, query_pt)

        expected_v = SVector(query_pt[2], query_pt[1], 0.0)
        @test v ≈ expected_v

        # J = [0 1 0; 1 0 0; 0 0 0]
        # Div = 0
        @test div ≈ 0.0 atol = 1e-10

        # Curl = (0,0,0) because J is symmetric
        @test norm(curl) ≈ 0.0 atol = 1e-10

        # S = 0.5(J + J') = J = [0 1 0; 1 0 0; 0 0 0]
        # σ = S - 0 = J
        expected_shear = SMatrix{3,3,Float64}(0, 1, 0, 1, 0, 0, 0, 0, 0)
        @test shear ≈ expected_shear
    end

    @testset "Bounds Checking" begin
        points_bad = make_grid_points(2) # 8 points
        velocities_bad = [SVector(0.0, 0.0, 0.0) for _ in 1:7] # Only 7 velocities

        @test_throws ArgumentError VelocityEstimator(points_bad, velocities_bad)

        est_den = DensityEstimator(points_bad)
        @test_throws ArgumentError VelocityEstimator(est_den, velocities_bad)
    end
end
