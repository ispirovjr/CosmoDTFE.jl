using Test
using JuliaDTFE
using StaticArrays

@testset "Composite Estimators" begin
    # 1. Basic check
    pointsBasic = [SVector{3,Float64}(rand(), rand(), rand()) for _ in 1:1000]
    weightsBasic = ones(Float64, 1000)

    @testset "Generation and basic evaluation" begin
        # use a tiny maxPoints to force branch creation
        compEst = CompositeEstimator(DensityEstimator, pointsBasic, weightsBasic, maxPoints=100, padding=0.1)

        # Test evaluation at a few points
        val = compEst(SVector(0.5, 0.5, 0.5))
        @test val >= 0.0 # Density is theoretically non-negative

        pt2 = SVector(0.0, 0.0, 0.0)
        val2 = compEst(pt2)
        @test val2 >= 0.0
    end

    @testset "VelocityEstimator generation" begin
        vels = [SVector{3,Float64}(1.0, 0.0, 0.0) for _ in 1:1000]
        compVel = CompositeEstimator(VelocityEstimator, pointsBasic, vels, maxPoints=100, padding=0.1)
        valVel = compVel(SVector(0.5, 0.5, 0.5))
        @test valVel isa SVector{3,Float64}

        vInterp, div, shear, vort = velocityGradient(compVel, SVector(0.5, 0.5, 0.5))
        @test vInterp isa SVector{3,Float64}
        @test div isa Float64
    end

    if get(ENV, "DO_HEAVY_TESTS", "false") == "true"
        @testset "Huge load survival" begin
            nPointsHeavy = 80_000_000
            @info "Running heavy test with $(nPointsHeavy) points..."
            pointsHeavy = [SVector{3,Float64}(rand(), rand(), rand()) for _ in 1:nPointsHeavy]
            weightsHeavy = ones(Float64, nPointsHeavy)

            # This shouldn't segfault
            compEstHeavy = CompositeEstimator(DensityEstimator, pointsHeavy, weightsHeavy, maxPoints=5_000_000, padding=0.01)
            @test compEstHeavy(SVector(0.5, 0.5, 0.5)) >= 0.0

            # too unstable  

            codeToRun = """
            using JuliaDTFE
            using StaticArrays
            pts = [SVector{3, Float64}(rand(), rand(), rand()) for _ in 1:$(nPointsHeavy)]
            w = ones(Float64, $(nPointsHeavy))
            est = DensityEstimator(pts, w) 
            """
            # We run it in a separate process so it doesn't kill the main test suite
            successBool = success(pipeline(`julia -e $codeToRun`))
            @test !successBool
        end
    else
        @info "Skipping heavy load tests. Set ENV[\"DO_HEAVY_TESTS\"]=\"true\" to run."
    end

    @testset "BVH Tree Leaf Count" begin
        Nside = 16
        pts = [SVector{3,Float64}(x, y, z) for x in 1:Nside for y in 1:Nside for z in 1:Nside]
        N = length(pts)
        Nmax = N ÷ 16
        weights = ones(Float64, N)
        compEst = CompositeEstimator(DensityEstimator, pts, weights, maxPoints=Nmax, padding=0.1)

        function countLeaves(node)
            if node isa CompositeBVHLeaf
                return 1
            else
                return countLeaves(node.leftChild) + countLeaves(node.rightChild)
            end
        end

        numLeaves = countLeaves(compEst.tree)
        @test numLeaves == 16
    end

    @testset "Parallel generation" begin
        Nside = 16
        pts = [SVector{3,Float64}(x, y, z) for x in 1:Nside for y in 1:Nside for z in 1:Nside]
        N = length(pts)
        Nmax = N ÷ 16
        weights = ones(Float64, N)

        function countLeaves(node)
            if node isa CompositeBVHLeaf
                return 1
            else
                return countLeaves(node.leftChild) + countLeaves(node.rightChild)
            end
        end

        # Build in parallel using all available threads
        compEstPar = CompositeEstimator(DensityEstimator, pts, weights, Threads.nthreads(), maxPoints=Nmax, padding=0.1)

        @test countLeaves(compEstPar.tree) == 16

        # Evaluation should return valid density
        val = compEstPar(SVector(8.0, 8.0, 8.0))
        @test val >= 0.0

        # Build serial for comparison
        compEstSer = CompositeEstimator(DensityEstimator, pts, weights, maxPoints=Nmax, padding=0.1)
        valSer = compEstSer(SVector(8.0, 8.0, 8.0))

        # Both should agree
        @test val ≈ valSer
    end

    @testset "Parallel VelocityEstimator" begin
        pts = [SVector{3,Float64}(rand(), rand(), rand()) for _ in 1:1000]
        vels = [SVector{3,Float64}(1.0, 0.0, 0.0) for _ in 1:1000]
        compVelPar = CompositeEstimator(VelocityEstimator, pts, vels, Threads.nthreads(), maxPoints=100, padding=0.1)
        valVel = compVelPar(SVector(0.5, 0.5, 0.5))
        @test valVel isa SVector{3,Float64}
    end

end

