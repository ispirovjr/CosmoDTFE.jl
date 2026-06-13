using Test
using CosmoDTFE
using LinearAlgebra

@testset "Composite Estimators" begin
    pointCloud(nSide=4) = [
        Point3(x, y, z)
        for x in range(0.0, 1.0, length=nSide)
        for y in range(0.0, 1.0, length=nSide)
        for z in range(0.0, 1.0, length=nSide)
    ]

    weightsFor(points) = [1.0 + 0.2 * point[1] + 0.1 * point[2] + 0.05 * point[3] for point in points]
    splitProbePoints = [
        Point3(0.48, 0.48, 0.48),
        Point3(0.52, 0.48, 0.48),
        Point3(0.48, 0.52, 0.48),
        Point3(0.48, 0.48, 0.52),
    ]

    function countLeaves(node)
        if node isa CompositeBVHLeaf
            return 1
        end
        return countLeaves(node.leftChild) + countLeaves(node.rightChild)
    end

    @testset "single-leaf density composite matches the global estimator" begin
        points = pointCloud()
        weights = weightsFor(points)
        globalEstimator = DensityEstimator(points, weights; depth=5)
        compositeEstimator = CompositeEstimator(DensityEstimator, points, weights; maxPoints=length(points), padding=0.0)

        queryPoints = [
            Point3(0.32, 0.35, 0.34),
            Point3(0.70, 0.45, 0.50),
            Point3(0.80, 0.82, 0.70),
        ]

        for queryPoint in queryPoints
            @test compositeEstimator(queryPoint) ≈ globalEstimator(queryPoint)
        end
    end

    @testset "split density composite is deterministic across serial and parallel builders" begin
        points = pointCloud()
        weights = weightsFor(points)
        maxPoints = div(length(points), 8)

        serialEstimator = CompositeEstimator(DensityEstimator, points, weights; maxPoints=maxPoints, padding=0.50)
        parallelEstimator = CompositeEstimator(DensityEstimator, points, weights, max(1, Threads.nthreads()); maxPoints=maxPoints, padding=0.50)

        @test countLeaves(serialEstimator.tree) == 8
        @test countLeaves(parallelEstimator.tree) == 8

        queryPoints = [
            Point3(0.35, 0.35, 0.35),
            Point3(0.80, 0.45, 0.55),
            Point3(0.88, 0.90, 0.85),
        ]

        for queryPoint in queryPoints
            @test parallelEstimator(queryPoint) ≈ serialEstimator(queryPoint)
        end
    end

    @testset "20 percent padded velocity composite matches global estimator near split planes" begin
        points = pointCloud(12)
        velocities = [Point3(point[2], -point[1], 0.25 * point[3]) for point in points]
        maxPoints = div(length(points), 8)

        globalEstimator = VelocityEstimator(points, velocities; depth=7)
        compositeEstimator = CompositeEstimator(VelocityEstimator, points, velocities; maxPoints=maxPoints, padding=0.20)

        for queryPoint in splitProbePoints
            @test compositeEstimator(queryPoint) ≈ globalEstimator(queryPoint)

            compVelocity, compDivergence, compShear, compVorticity = velocityGradient(compositeEstimator, queryPoint)
            globalVelocity, globalDivergence, globalShear, globalVorticity = velocityGradient(globalEstimator, queryPoint)

            @test compVelocity ≈ globalVelocity
            @test compDivergence ≈ globalDivergence
            @test compShear ≈ globalShear
            @test compVorticity ≈ globalVorticity
        end
    end

    @testset "20 percent padded density composite matches global estimator near split planes on a resolved grid" begin
        points = pointCloud(24)
        weights = weightsFor(points)
        maxPoints = div(length(points), 8)

        globalEstimator = DensityEstimator(points, weights; depth=4)
        compositeEstimator = CompositeEstimator(DensityEstimator, points, weights; maxPoints=maxPoints, padding=0.20)

        @test countLeaves(compositeEstimator.tree) == 8

        for queryPoint in splitProbePoints
            @test compositeEstimator(queryPoint) ≈ globalEstimator(queryPoint) rtol = 1e-10 atol = 1e-10
        end
    end

    @testset "constant velocity field survives composite routing" begin
        points = pointCloud()
        constantVelocity = Point3(1.0, -2.0, 0.5)
        velocities = [constantVelocity for _ in points]
        estimator = CompositeEstimator(VelocityEstimator, points, velocities; maxPoints=div(length(points), 8), padding=0.50)

        for queryPoint in (points[10], points[32], points[55])
            velocityValue = estimator(queryPoint)
            velocityInterp, divergence, shear, vorticity = velocityGradient(estimator, queryPoint)

            @test velocityValue ≈ constantVelocity
            @test velocityInterp ≈ constantVelocity
            @test divergence ≈ 0.0 atol = 1e-10
            @test norm(shear) ≈ 0.0 atol = 1e-10
            @test norm(vorticity) ≈ 0.0 atol = 1e-10
        end
    end

    @testset "parallel velocity composite matches serial composite" begin
        points = pointCloud()
        velocities = [Point3(point[2], -point[1], 0.25 * point[3]) for point in points]
        maxPoints = div(length(points), 8)

        serialEstimator = CompositeEstimator(VelocityEstimator, points, velocities; maxPoints=maxPoints, padding=0.50)
        parallelEstimator = CompositeEstimator(VelocityEstimator, points, velocities, max(1, Threads.nthreads()); maxPoints=maxPoints, padding=0.50)

        for queryPoint in (points[10], points[32], points[55])
            @test parallelEstimator(queryPoint) ≈ serialEstimator(queryPoint)
            @test velocityGradient(parallelEstimator, queryPoint)[1] ≈ velocityGradient(serialEstimator, queryPoint)[1]
        end
    end
end
