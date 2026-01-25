using Test
using StaticArrays
using JuliaDTFE

@testset "Estimators" begin

    @testset "standardEstimator builds valid structure" begin
        # Create simple point cloud
        points = [Point3(rand(), rand(), rand()) for _ in 1:50]

        bvh, triangulation, tets = standardEstimator(points, 6)

        @test bvh isa BoundingVolumeHierarchy
        @test triangulation isa Triangulation3D
        @test length(triangulation.points) == 50
        @test length(triangulation.rhoStar) == 50
        @test size(tets, 2) == 4  # Each tet has 4 vertices
    end

    @testset "standardEstimator with weights" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        weights = rand(30)

        bvh, triangulation, tets = standardEstimator(points, weights, 5)

        @test bvh isa BoundingVolumeHierarchy
        @test length(triangulation.rhoStar) == 30
        # Densities should vary with weights
        @test !all(triangulation.rhoStar .== triangulation.rhoStar[1])
    end

    @testset "dtfe returns density at point" begin
        # Create uniform point cloud in unit cube
        n = 5
        points = [Point3(x, y, z) for x in 0:0.25:1 for y in 0:0.25:1 for z in 0:0.25:1]

        bvh, triangulation, tets = standardEstimator(points, 6)

        # Test point inside domain
        testPoint = [0.5, 0.5, 0.5]
        density = dtfe(testPoint, bvh, tets, triangulation)

        @test density >= 0.0
        @test isfinite(density)
    end

    @testset "dtfe returns 0 outside domain" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        bvh, triangulation, tets = standardEstimator(points, 5)

        # Point far outside
        outsidePoint = [100.0, 100.0, 100.0]
        density = dtfe(outsidePoint, bvh, tets, triangulation)

        @test density == 0.0
    end

    @testset "dtfeMultiThread produces 3D grid" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        bvh, triangulation, tets = standardEstimator(points, 5)

        # Small grid for testing
        xs = range(0.0, 1.0, length=3)
        ys = range(0.0, 1.0, length=3)
        zs = range(0.0, 1.0, length=3)

        densityGrid = dtfeMultiThread((xs, ys, zs), bvh, tets, triangulation)

        @test size(densityGrid) == (3, 3, 3)
        @test all(densityGrid .>= 0.0)
    end

end
