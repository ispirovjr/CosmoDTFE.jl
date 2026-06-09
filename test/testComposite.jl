using Test
using CosmoDTFE
using StaticArrays

@testset "Composite Estimators" begin
    pointsBasic = [SVector{3,Float64}(rand(), rand(), rand()) for _ in 1:1000]
    weightsBasic = ones(Float64, 1000)

    @testset "Generation and basic evaluation" begin
        compEst = CompositeEstimator(DensityEstimator, pointsBasic, weightsBasic, maxPoints=100, padding=0.1)

        valueOne = compEst(SVector(0.5, 0.5, 0.5))
        @test valueOne >= 0.0

        valueTwo = compEst(SVector(0.0, 0.0, 0.0))
        @test valueTwo >= 0.0
    end

    @testset "VelocityEstimator generation" begin
        velocities = [SVector{3,Float64}(1.0, 0.0, 0.0) for _ in 1:1000]
        compVel = CompositeEstimator(VelocityEstimator, pointsBasic, velocities, maxPoints=100, padding=0.1)
        velocityValue = compVel(SVector(0.5, 0.5, 0.5))
        @test velocityValue isa SVector{3,Float64}

        velocityInterp, divergence, shear, vorticity = velocityGradient(compVel, SVector(0.5, 0.5, 0.5))
        @test velocityInterp isa SVector{3,Float64}
        @test divergence isa Float64
        @test shear isa SMatrix{3,3,Float64}
        @test vorticity isa SVector{3,Float64}
    end

    if get(ENV, "DO_HEAVY_TESTS", "false") == "true"
        @testset "Huge load survival" begin
            nPointsHeavy = 80_000_000
            @info "Running heavy test with $(nPointsHeavy) points..."
            pointsHeavy = [SVector{3,Float64}(rand(), rand(), rand()) for _ in 1:nPointsHeavy]
            weightsHeavy = ones(Float64, nPointsHeavy)

            compEstHeavy = CompositeEstimator(DensityEstimator, pointsHeavy, weightsHeavy, maxPoints=5_000_000, padding=0.01)
            @test compEstHeavy(SVector(0.5, 0.5, 0.5)) >= 0.0

            codeToRun = """
            using CosmoDTFE
            using StaticArrays
            points = [SVector{3, Float64}(rand(), rand(), rand()) for _ in 1:$(nPointsHeavy)]
            weights = ones(Float64, $(nPointsHeavy))
            estimator = DensityEstimator(points, weights)
            """
            successBool = success(pipeline(`julia -e $codeToRun`))
            @test !successBool
        end
    else
        @info "Skipping heavy load tests. Set ENV[\"DO_HEAVY_TESTS\"]=\"true\" to run."
    end

    @testset "BVH Tree Leaf Count" begin
        nSide = 16
        points = [SVector{3,Float64}(x, y, z) for x in 1:nSide for y in 1:nSide for z in 1:nSide]
        pointCount = length(points)
        maxPoints = div(pointCount, 16)
        weights = ones(Float64, pointCount)
        compEst = CompositeEstimator(DensityEstimator, points, weights, maxPoints=maxPoints, padding=0.1)

        function countLeaves(node)
            if node isa CompositeBVHLeaf
                return 1
            end
            return countLeaves(node.leftChild) + countLeaves(node.rightChild)
        end

        @test countLeaves(compEst.tree) == 16
    end

    @testset "Parallel generation" begin
        nSide = 16
        points = [SVector{3,Float64}(x, y, z) for x in 1:nSide for y in 1:nSide for z in 1:nSide]
        pointCount = length(points)
        maxPoints = div(pointCount, 16)
        weights = ones(Float64, pointCount)

        function countLeaves(node)
            if node isa CompositeBVHLeaf
                return 1
            end
            return countLeaves(node.leftChild) + countLeaves(node.rightChild)
        end

        compEstPar = CompositeEstimator(DensityEstimator, points, weights, Threads.nthreads(), maxPoints=maxPoints, padding=0.1)
        @test countLeaves(compEstPar.tree) == 16

        valuePar = compEstPar(SVector(8.0, 8.0, 8.0))
        @test valuePar >= 0.0

        compEstSer = CompositeEstimator(DensityEstimator, points, weights, maxPoints=maxPoints, padding=0.1)
        valueSer = compEstSer(SVector(8.0, 8.0, 8.0))
        @test valuePar ≈ valueSer
    end

    @testset "Parallel VelocityEstimator" begin
        points = [SVector{3,Float64}(rand(), rand(), rand()) for _ in 1:1000]
        velocities = [SVector{3,Float64}(1.0, 0.0, 0.0) for _ in 1:1000]
        compVelPar = CompositeEstimator(VelocityEstimator, points, velocities, Threads.nthreads(), maxPoints=100, padding=0.1)
        velocityValue = compVelPar(SVector(0.5, 0.5, 0.5))
        @test velocityValue isa SVector{3,Float64}
    end
end
