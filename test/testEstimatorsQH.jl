using Test
using StaticArrays
using JuliaDTFE
using Statistics
using LinearAlgebra

# Helper: build a DensityEstimator using the Quickhull backend
function makeDensityEstimatorQH(points, depth::Int=9)
    coords, tets = tessellateQH(points)
    triangulation = Triangulation3D(points, tets)

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    return DensityEstimator(bvh, triangulation, tets)
end

function makeDensityEstimatorQH(points, weights::Vector, depth::Int=9)
    coords, tets = tessellateQH(points)
    triangulation = Triangulation3D(points, tets, weights)

    simplices = coords[:, tets]
    bvh = BoundingVolumeHierarchy(simplices, depth)

    return DensityEstimator(bvh, triangulation, tets)
end


@testset "EstimatorsQH" begin

    @testset "DensityEstimator (QH) builds valid structure" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:50]

        estimator = makeDensityEstimatorQH(points, 6)

        @test estimator isa DensityEstimator
        @test estimator.bvh isa BoundingVolumeHierarchy
        @test estimator.triangulation isa Triangulation3D
        @test length(estimator.triangulation.points) == 50
        @test length(estimator.triangulation.rhoStar) == 50
        @test size(estimator.tetrahedra, 2) == 4
    end

    @testset "DensityEstimator (QH) with weights" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        weights = rand(30)

        estimator = makeDensityEstimatorQH(points, weights, 5)

        @test estimator isa DensityEstimator
        @test length(estimator.triangulation.rhoStar) == 30
        @test !all(estimator.triangulation.rhoStar .== estimator.triangulation.rhoStar[1])
    end

    @testset "dtfe (QH) returns density at point" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:100]

        estimator = makeDensityEstimatorQH(points, 6)

        testPoint = [0.5, 0.5, 0.5]
        density = estimator(testPoint)

        @test density >= 0.0
        @test isfinite(density)
    end

    @testset "dtfe (QH) with deterministic spherical grid - analytical check" begin

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

        estimator = makeDensityEstimatorQH(points, 6)

        # Mass Conservation
        totalMass = sum(begin
            pos = estimator.triangulation.points[tet]
            vol = JuliaDTFE.computeVolume(pos)
            rhos = estimator.triangulation.rhoStar[tet]
            mean(rhos) * vol
        end for tet in eachrow(estimator.tetrahedra))

        @test totalMass ≈ N

        # Density Scaling
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

    @testset "dtfe (QH) returns 0 outside domain" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        estimator = makeDensityEstimatorQH(points, 5)

        outsidePoint = [100.0, 100.0, 100.0]
        density = estimator(outsidePoint)

        @test density == 0.0
    end

    @testset "dtfeMultiThread (QH) produces 3D grid" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        estimator = makeDensityEstimatorQH(points, 5)

        xs = range(0.0, 1.0, length=3)
        ys = range(0.0, 1.0, length=3)
        zs = range(0.0, 1.0, length=3)

        densityGrid = estimator((xs, ys, zs))

        @test size(densityGrid) == (3, 3, 3)
        @test all(densityGrid .>= 0.0)
    end

    @testset "dtfeMultiThread (QH) produces 1D array" begin
        points = [Point3(rand(), rand(), rand()) for _ in 1:30]
        estimator = makeDensityEstimatorQH(points, 5)

        query_points = [Point3(rand(), rand(), rand()) for _ in 1:10]
        densityArray = estimator(query_points)

        @test length(densityArray) == 10
        @test all(densityArray .>= 0.0)
    end

end
