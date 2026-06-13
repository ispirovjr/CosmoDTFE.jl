using Test
using LinearAlgebra
using CosmoDTFE

@testset "PeriodicEstimator" begin
    gridPoints(boxMin=Point3(0.0, 0.0, 0.0), boxSize=Point3(1.0, 1.0, 1.0), nSide=4) = [
        boxMin + Point3(fx * boxSize[1], fy * boxSize[2], fz * boxSize[3])
        for fx in range(0.05, 0.95, length=nSide)
        for fy in range(0.05, 0.95, length=nSide)
        for fz in range(0.05, 0.95, length=nSide)
    ]

    @testset "periodicCopies duplicates faces, edges, and corners exactly" begin
        points = [
            Point3(0.50, 0.50, 0.50),
            Point3(0.05, 0.50, 0.50),
            Point3(0.05, 0.05, 0.50),
            Point3(0.05, 0.05, 0.05),
            Point3(0.95, 0.95, 0.95),
        ]
        ids = collect(1:length(points))

        copyPoints, copyIds = CosmoDTFE.periodicCopies(points, ids, Point3(0.0, 0.0, 0.0), Point3(1.0, 1.0, 1.0), 0.10)

        @test length(copyPoints) == 23
        @test count(==(1), copyIds) == 1
        @test count(==(2), copyIds) == 2
        @test count(==(3), copyIds) == 4
        @test count(==(4), copyIds) == 8
        @test count(==(5), copyIds) == 8
        @test any(point -> isapprox(point, Point3(1.05, 0.50, 0.50); atol=1e-14, rtol=0.0), copyPoints)
        @test any(point -> isapprox(point, Point3(-0.05, -0.05, -0.05); atol=1e-14, rtol=0.0), copyPoints)
    end

    @testset "periodicCopies validates property length" begin
        points = [Point3(0.05, 0.50, 0.50), Point3(0.50, 0.50, 0.50)]

        @test_throws ArgumentError CosmoDTFE.periodicCopies(points, [1], Point3(0.0, 0.0, 0.0), Point3(1.0, 1.0, 1.0), 0.10)
    end

    @testset "density wrapping is invariant across all box directions" begin
        boxMin = Point3(-1.0, 2.0, 4.0)
        boxSize = Point3(2.0, 3.0, 5.0)
        points = gridPoints(boxMin, boxSize, 4)
        weights = ones(Float64, length(points))
        estimator = PeriodicEstimator(DensityEstimator, points, weights; boxMin=boxMin, boxSize=boxSize, padding=0.12, depth=5)

        queryPoint = boxMin + Point3(0.02 * boxSize[1], 0.50 * boxSize[2], 0.50 * boxSize[3])
        densityValue = estimator(queryPoint)

        @test densityValue > 0.0
        @test isfinite(densityValue)
        @test estimator(queryPoint + Point3(boxSize[1], 0.0, 0.0)) ≈ densityValue
        @test estimator(queryPoint - Point3(0.0, boxSize[2], 0.0)) ≈ densityValue
        @test estimator(queryPoint + Point3(0.0, 0.0, 2.0 * boxSize[3])) ≈ densityValue
        @test estimator(queryPoint + boxSize) ≈ densityValue
    end

    @testset "unweighted density constructor is periodic" begin
        points = gridPoints()
        estimator = PeriodicEstimator(DensityEstimator, points; boxSize=Point3(1.0, 1.0, 1.0), padding=0.12, depth=5)

        @test estimator(Point3(-0.02, 0.50, 0.50)) ≈ estimator(Point3(0.98, 0.50, 0.50))
    end

    @testset "velocity wrapping preserves a constant field and zero gradient" begin
        boxMin = Point3(-1.0, 2.0, 4.0)
        boxSize = Point3(2.0, 3.0, 5.0)
        points = gridPoints(boxMin, boxSize, 4)
        constantVelocity = Point3(2.0, -1.0, 0.5)
        velocities = [constantVelocity for _ in points]
        estimator = PeriodicEstimator(VelocityEstimator, points, velocities; boxMin=boxMin, boxSize=boxSize, padding=0.12, depth=5)

        queryPoint = boxMin + Point3(1.02 * boxSize[1], 0.50 * boxSize[2], 0.50 * boxSize[3])
        velocityValue = estimator(queryPoint)
        velocityInterp, divergence, shear, vorticity = velocityGradient(estimator, queryPoint - Point3(boxSize[1], 0.0, 0.0))

        @test velocityValue ≈ constantVelocity
        @test velocityInterp ≈ constantVelocity
        @test divergence ≈ 0.0 atol = 1e-10
        @test norm(shear) ≈ 0.0 atol = 1e-10
        @test norm(vorticity) ≈ 0.0 atol = 1e-10
    end
end
