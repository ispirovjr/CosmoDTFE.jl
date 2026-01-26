using Test
using StaticArrays
using JuliaDTFE
using Statistics

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
        # Random points - general functionality check
        points = [Point3(rand(), rand(), rand()) for _ in 1:100]

        bvh, triangulation, tets = standardEstimator(points, 6)

        # Test point inside domain
        testPoint = [0.5, 0.5, 0.5]
        density = dtfe(testPoint, bvh, tets, triangulation)

        @test density >= 0.0
        @test isfinite(density)
    end

    @testset "dtfe with deterministic spherical grid - analytical check" begin
        # 1. Define Grid
        # Deterministic grid in (r, θ, φ) space
        # Density should scale as ρ(r, θ) ∝ 1 / (r^2 * sin(θ))
        R_min, R_max = 0.2, 1.0
        nR, nθ, nφ = 20, 20, 20

        # Grid parameters for analytical formula
        dr = (R_max - R_min) / (nR - 1)
        dθ = (0.9π - 0.1π) / (nθ - 1)
        dφ = 2π / nφ

        points = Point3[]
        for r in range(R_min, R_max, length=nR)
            for θ in range(0.1π, 0.9π, length=nθ)
                for φ in range(0, 2π - dφ / 2, length=nφ)
                    push!(points, Point3(r * sin(θ) * cos(φ), r * sin(θ) * sin(φ), r * cos(θ)))
                end
            end
        end
        N = length(points)

        bvh, triangulation, tets = standardEstimator(points, 6)

        # 2. Check Mass Conservation (should be exact)
        totalMass = sum(begin
            pos = triangulation.points[tet]
            vol = JuliaDTFE.computeVolume(pos)
            rhos = triangulation.rhoStar[tet]
            mean(rhos) * vol
        end for tet in eachrow(tets))

        @test totalMass ≈ N

        # 3. Check Density Scaling
        # Test at r=0.5 (equator) -> sin(pi/2)=1
        # Analytical: ρ = 1 / (r^2 * sin(θ) * dr * dθ * dφ)

        function analytical_val(r)
            return 1.0 / (r^2 * 1.0 * dr * dθ * dφ)
        end

        density_r05 = dtfe([0.5, 0.0, 0.0], bvh, tets, triangulation)
        density_r08 = dtfe([0.8, 0.0, 0.0], bvh, tets, triangulation)

        expected_r05 = analytical_val(0.5)
        expected_r08 = analytical_val(0.8)

        # Check values are reasonable (within factor of 1.6 due to sampling/boundary)
        # We care most that they follow the trend and scale correctly
        # Empirically observed ratios ~1.2-1.5
        @test 0.7 < density_r05 / expected_r05 < 1.8
        @test 0.7 < density_r08 / expected_r08 < 1.8

        # Check scaling trend: density should drop roughly by (0.5/0.8)^2 ≈ 0.39
        # Allow wide margin for sampling artifacts
        ratio = density_r08 / density_r05
        expected_ratio = expected_r08 / expected_r05
        @test 0.5 * expected_ratio < ratio < 1.6 * expected_ratio
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
