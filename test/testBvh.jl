using Test
using StaticArrays
using JuliaDTFE

@testset "BVH" begin

    @testset "BVHLeaf construction" begin
        # Test leaf node construction
        leaf = BVHLeaf([1, 2, 3])
        @test leaf.data == [1, 2, 3]
        @test leaf isa BVHTree
    end

    @testset "BVHNode construction" begin
        # Test node construction with children
        leftLeaf = BVHLeaf([1, 2])
        rightLeaf = BVHLeaf([3, 4])
        node = BVHNode(3, leftLeaf, rightLeaf)

        @test node.depth == 3
        @test node.leftChild === leftLeaf
        @test node.rightChild === rightLeaf
        @test node isa BVHTree
    end

    @testset "BoundingVolumeHierarchy from array" begin
        # Create simple test simplices (2 tetrahedra as 3D array)
        # Shape: [3 coords, 4 vertices, N simplices]
        simplices = zeros(Float64, 3, 4, 2)

        # First tetrahedron at origin
        simplices[:, 1, 1] = [0.0, 0.0, 0.0]
        simplices[:, 2, 1] = [1.0, 0.0, 0.0]
        simplices[:, 3, 1] = [0.0, 1.0, 0.0]
        simplices[:, 4, 1] = [0.0, 0.0, 1.0]

        # Second tetrahedron offset
        simplices[:, 1, 2] = [1.0, 1.0, 1.0]
        simplices[:, 2, 2] = [2.0, 1.0, 1.0]
        simplices[:, 3, 2] = [1.0, 2.0, 1.0]
        simplices[:, 4, 2] = [1.0, 1.0, 2.0]

        bvh = BoundingVolumeHierarchy(simplices, 3)

        @test bvh isa BoundingVolumeHierarchy
        @test bvh.tree isa BVHTree
        @test size(bvh.bbox) == (3, 2)  # 3D bounding box
    end

    @testset "BoundingVolumeHierarchy depth handling" begin
        # Test that depth=0 creates leaf nodes
        simplices = zeros(Float64, 3, 4, 1)
        simplices[:, 1, 1] = [0.0, 0.0, 0.0]
        simplices[:, 2, 1] = [1.0, 0.0, 0.0]
        simplices[:, 3, 1] = [0.0, 1.0, 0.0]
        simplices[:, 4, 1] = [0.0, 0.0, 1.0]

        bvh = BoundingVolumeHierarchy(simplices, 0)
        @test bvh.tree isa BVHLeaf
    end

end
