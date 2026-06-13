using Test
using CosmoDTFE

@testset "BVH" begin
    simplices = zeros(Float64, 3, 3, 4)

    simplices[:, 1, 1] = [0.0, 0.0, 0.0]
    simplices[:, 1, 2] = [1.0, 0.0, 0.0]
    simplices[:, 1, 3] = [0.0, 1.0, 0.0]
    simplices[:, 1, 4] = [0.0, 0.0, 1.0]

    simplices[:, 2, 1] = [2.0, 0.0, 0.0]
    simplices[:, 2, 2] = [3.0, 0.0, 0.0]
    simplices[:, 2, 3] = [2.0, 1.0, 0.0]
    simplices[:, 2, 4] = [2.0, 0.0, 1.0]

    simplices[:, 3, 1] = [0.25, 0.25, 0.25]
    simplices[:, 3, 2] = [0.75, 0.25, 0.25]
    simplices[:, 3, 3] = [0.25, 0.75, 0.25]
    simplices[:, 3, 4] = [0.25, 0.25, 0.75]

    @testset "global bounding box matches simplex extents" begin
        bvh = BoundingVolumeHierarchy(simplices, 3)

        @test bvh.bbox[:, 1] ≈ [0.0, 0.0, 0.0]
        @test bvh.bbox[:, 2] ≈ [3.0, 1.0, 1.0]
    end

    @testset "depth zero stores all simplex ids in one leaf" begin
        bvh = BoundingVolumeHierarchy(simplices, 0)

        @test bvh.tree isa BVHLeaf
        @test bvh.tree.data == [1, 2, 3]
    end

    @testset "split-spanning boxes stay reachable from both sides" begin
        simplices = zeros(Float64, 3, 1, 4)
        simplices[:, 1, 1] = [0.25, 0.25, 0.25]
        simplices[:, 1, 2] = [0.75, 0.25, 0.25]
        simplices[:, 1, 3] = [0.25, 0.75, 0.25]
        simplices[:, 1, 4] = [0.25, 0.25, 0.75]

        box = [0.0 1.0; 0.0 1.0; 0.0 1.0]
        bvh = BoundingVolumeHierarchy(simplices, 1; box=box)

        @test 1 in recursiveSearch(Point3(0.25, 0.25, 0.25), bvh.tree, bvh.bbox)
        @test 1 in recursiveSearch(Point3(0.25, 0.75, 0.25), bvh.tree, bvh.bbox)
    end

    @testset "leaf equality ignores candidate order" begin
        @test BVHLeaf([1, 2, 3]) == BVHLeaf([3, 2, 1])
        @test BVHLeaf([1, 2, 3]) != BVHLeaf([1, 2])
    end
end
