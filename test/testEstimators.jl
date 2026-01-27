using Test
using StaticArrays
using JuliaDTFE
using Statistics

@testset "Estimators" begin

    @testset "DensityEstimator builds valid structure" begin
        # Create simple point cloud
        points = [Point3(rand(), rand(), rand()) for _ in 1:50]

        estimator = DensityEstimator(points, 6)

        @test estimator isa DensityEstimator
        @test estimator.bvh isa BoundingVolumeHierarchy
        @test estimator.triangulation isa Triangulation3D
        @test length(estimator.triangulation.points) == 50
        @test length(estimator.triangulation.rhoStar) == 50
        @test size(estimator.tetrahedra, 2) == 4  # Each tet has 4 vertices
    end

    @testset "DensityEstimator with weights" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        weights = rand(30)

        estimator = DensityEstimator(points, weights, 5)

        @test estimator isa DensityEstimator
        @test length(estimator.triangulation.rhoStar) == 30
        # Densities should vary with weights
        @test !all(estimator.triangulation.rhoStar .== estimator.triangulation.rhoStar[1])
    end

    @testset "dtfe returns density at point" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:100]

        estimator = DensityEstimator(points, 6)

        testPoint = [0.5, 0.5, 0.5]
        # Use functor call
        density = estimator(testPoint)

        @test density >= 0.0
        @test isfinite(density)
    end

    @testset "dtfe with deterministic spherical grid - analytical check" begin

        # Density should scale as ρ(r, θ) ∝ 1 / (r^2 * sin(θ))
        rMin, rMax = 0.2, 1.0
        nR, nθ, nφ = 30, 30, 30

        dr = (rMax - rMin) / (nR - 1)
        dθ = (0.9π - 0.1π) / (nθ - 1)
        dφ = 2π / nφ

        points = Point3[]
        for r in range(rMin, rMax, length=nR)
            for θ in range(0.1π, 0.9π, length=nθ)
                for φ in range(0, 2π - dφ / 2, length=nφ)
                    push!(points, Point3(r * sin(θ) * cos(φ), r * sin(θ) * sin(φ), r * cos(θ)))
                end
            end
        end
        N = length(points)

        estimator = DensityEstimator(points, 6)

        # 2. Check Mass Conservation
        totalMass = sum(begin
            pos = estimator.triangulation.points[tet]
            vol = JuliaDTFE.computeVolume(pos)
            rhos = estimator.triangulation.rhoStar[tet]
            mean(rhos) * vol
        end for tet in eachrow(estimator.tetrahedra))

        @test totalMass ≈ N

        # 3. Check Density Scaling
        function analyticalVal(r)
            return 1.0 / (r^2 * 1.0 * dr * dθ * dφ)
        end

        densityR05 = estimator([0.5, 0.0, 0.0])
        densityR08 = estimator([0.8, 0.0, 0.0])

        expectedR05 = analyticalVal(0.5)
        expectedR08 = analyticalVal(0.8)

        @test 0.7 < densityR05 / expectedR05 < 1.8
        @test 0.7 < densityR08 / expectedR08 < 1.8

    end

    @testset "dtfe returns 0 outside domain" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        estimator = DensityEstimator(points, 5)

        # Point far outside
        outsidePoint = [100.0, 100.0, 100.0]
        density = estimator(outsidePoint)

        @test density == 0.0
    end

    @testset "dtfeMultiThread produces 3D grid" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        estimator = DensityEstimator(points, 5)

        # Small grid for testing
        xs = range(0.0, 1.0, length=3)
        ys = range(0.0, 1.0, length=3)
        zs = range(0.0, 1.0, length=3)

        # Use functor call on tuple
        densityGrid = estimator((xs, ys, zs))

        @test size(densityGrid) == (3, 3, 3)
        @test all(densityGrid .>= 0.0)
    end

    @testset "dtfeMultiThread produces 1D array" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        estimator = DensityEstimator(points, 5)

        # List of points for testing
        query_points = [Point3(rand(), rand(), rand()) for _ in 1:10]

        # Use functor call on vector
        densityArray = estimator(query_points)

        @test length(densityArray) == 10
        @test all(densityArray .>= 0.0)
    end

end
