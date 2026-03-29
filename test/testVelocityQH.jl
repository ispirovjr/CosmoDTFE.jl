using Test
using StaticArrays
using LinearAlgebra
using JuliaDTFE

# Helper: build a VelocityEstimator using the Quickhull backend
function makeVelocityEstimatorQH(points, velocities::Vector, depth::Int=9)
    coords, tets = tessellateQH(points)
    triangulation = Triangulation3D(points, tets)

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    vels = [SVector{3,Float64}(v) for v in velocities]
    return VelocityEstimator(bvh, triangulation, tets, vels)
end

# Setup a simple 3D grid of points to tessellate
function makeGridPointsQH(N=5)
    points = [SVector{3,Float64}(x, y, z) for x in 1:N, y in 1:N, z in 1:N]
    return vec(points)
end

@testset "VelocityEstimatorQH Tests" begin
    points = makeGridPointsQH(5) # 5x5x5 grid

    @testset "Constant Velocity Field (QH)" begin
        const_vel = SVector(1.0, 2.0, 3.0)
        velocities = [const_vel for _ in points]

        est = makeVelocityEstimatorQH(points, velocities)

        query_pt = SVector(2.5, 2.5, 2.5)

        v, div, shear, curl = velocityGradient(est, query_pt)

        @test v ≈ const_vel
        @test div ≈ 0.0 atol = 1e-10
        @test norm(shear) ≈ 0.0 atol = 1e-10
        @test norm(curl) ≈ 0.0 atol = 1e-10
    end

    @testset "Pure Expansion (Linear Field) (QH)" begin
        velocities = [p for p in points]

        est = makeVelocityEstimatorQH(points, velocities)
        query_pt = SVector(2.5, 2.5, 2.5)

        v, div, shear, curl = velocityGradient(est, query_pt)

        @test v ≈ query_pt
        @test div ≈ 3.0
        @test norm(shear) ≈ 0.0 atol = 1e-10
        @test norm(curl) ≈ 0.0 atol = 1e-10
    end

    @testset "Solid Body Rotation (QH)" begin
        velocities = [SVector(-p[2], p[1], 0.0) for p in points]

        est = makeVelocityEstimatorQH(points, velocities)
        query_pt = SVector(2.5, 2.5, 2.5)

        v, div, shear, curl = velocityGradient(est, query_pt)

        expected_v = SVector(-query_pt[2], query_pt[1], 0.0)
        @test v ≈ expected_v
        @test div ≈ 0.0 atol = 1e-10
        @test curl ≈ SVector(0.0, 0.0, 2.0)
        @test norm(shear) ≈ 0.0 atol = 1e-10
    end

    @testset "Pure Shear (QH)" begin
        velocities = [SVector(p[2], p[1], 0.0) for p in points]

        est = makeVelocityEstimatorQH(points, velocities)
        query_pt = SVector(2.5, 2.5, 2.5)

        v, div, shear, curl = velocityGradient(est, query_pt)

        expected_v = SVector(query_pt[2], query_pt[1], 0.0)
        @test v ≈ expected_v
        @test div ≈ 0.0 atol = 1e-10
        @test norm(curl) ≈ 0.0 atol = 1e-10

        expected_shear = SMatrix{3,3,Float64}(0, 1, 0, 1, 0, 0, 0, 0, 0)
        @test shear ≈ expected_shear
    end
end
