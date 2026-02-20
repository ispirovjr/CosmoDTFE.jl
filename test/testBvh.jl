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

        simplices = zeros(Float64, 3, 2, 4)

        # First tetrahedron at origin
        simplices[:, 1, 1] = [0.0, 0.0, 0.0]
        simplices[:, 1, 2] = [1.0, 0.0, 0.0]
        simplices[:, 1, 3] = [0.0, 1.0, 0.0]
        simplices[:, 1, 4] = [0.0, 0.0, 1.0]

        # Second tetrahedron offset
        simplices[:, 2, 1] = [1.0, 1.0, 1.0]
        simplices[:, 2, 2] = [2.0, 1.0, 1.0]
        simplices[:, 2, 3] = [1.0, 2.0, 1.0]
        simplices[:, 2, 4] = [1.0, 1.0, 2.0]

        bvh = BoundingVolumeHierarchy(simplices, 3)

        @test bvh isa BoundingVolumeHierarchy
        @test bvh.tree isa BVHTree
        @test size(bvh.bbox) == (3, 2)  # 3D bounding box

    end

    @testset "BoundingVolumeHierarchy depth handling" begin
        # Test that depth=0 creates leaf nodes
        simplices = zeros(Float64, 3, 1, 4)
        simplices[:, 1, 1] = [0.0, 0.0, 0.0]
        simplices[:, 1, 2] = [1.0, 0.0, 0.0]
        simplices[:, 1, 3] = [0.0, 1.0, 0.0]
        simplices[:, 1, 4] = [0.0, 0.0, 1.0]

        bvh = BoundingVolumeHierarchy(simplices, 0)
        @test bvh.tree isa BVHLeaf
    end


    @testset "Equality Operators" begin
        # Leaf equality
        l1 = BVHLeaf([1, 2, 3])
        l2 = BVHLeaf([3, 2, 1])
        l3 = BVHLeaf([1, 2])
        @test l1 == l2
        @test l1 != l3

        # Node equality
        n1 = BVHNode(1, l1, l3)
        n2 = BVHNode(1, l2, l3) # Identical structure since l1==l2
        n3 = BVHNode(1, l3, l1) # Different structure
        @test n1 == n2
        @test n1 != n3
    end

    @testset "Leaf Uniqueness Verification" begin
        # Test if given a sufficiently deep BVH, all leaves are unique (<= 1 element)

        # Create 3 tetrahedra
        # Create 3 tetrahedra
        simplices = zeros(Float64, 3, 3, 4)
        simplices[:, 1, :] .= [0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0; 0.0 1.0 0.0 1.0] # Reshaped
        # Wait, the broadcasting is safer with simple assignment
        simplices[:, 1, 1] = [0.0, 0.0, 0.0]
        simplices[:, 1, 2] = [1.0, 0.0, 0.0]
        simplices[:, 1, 3] = [0.0, 1.0, 0.0]
        simplices[:, 1, 4] = [0.0, 0.0, 1.0]
        simplices[:, 2, 1] = [2.0, 0.0, 0.0]
        simplices[:, 2, 2] = [3.0, 0.0, 0.0]
        simplices[:, 2, 3] = [2.0, 1.0, 0.0]
        simplices[:, 2, 4] = [2.0, 0.0, 1.0]
        simplices[:, 3, 1] = [4.0, 0.0, 0.0]
        simplices[:, 3, 2] = [5.0, 0.0, 0.0]
        simplices[:, 3, 3] = [4.0, 1.0, 0.0]
        simplices[:, 3, 4] = [4.0, 0.0, 1.0]

        bvh = BoundingVolumeHierarchy(simplices, 6)

        function verifyLeafProperties(node::BVHLeaf)
            @test length(node.data) <= 1
        end
        function verifyLeafProperties(node::BVHNode)
            verifyLeafProperties(node.leftChild)
            verifyLeafProperties(node.rightChild)
        end

        verifyLeafProperties(bvh.tree)
    end

    @testset "Node Structure Verification" begin

        simplices = zeros(Float64, 3, 3, 4)
        simplices[:, 1, 1] = [0.0, 0.0, 0.0]
        simplices[:, 1, 2] = [1.0, 0.0, 0.0]
        simplices[:, 1, 3] = [0.0, 1.0, 0.0]
        simplices[:, 1, 4] = [0.0, 0.0, 1.0]
        simplices[:, 2, 1] = [2.0, 0.0, 0.0]
        simplices[:, 2, 2] = [3.0, 0.0, 0.0]
        simplices[:, 2, 3] = [2.0, 1.0, 0.0]
        simplices[:, 2, 4] = [2.0, 0.0, 1.0]
        simplices[:, 3, 1] = [4.0, 0.0, 0.0]
        simplices[:, 3, 2] = [5.0, 0.0, 0.0]
        simplices[:, 3, 3] = [4.0, 1.0, 0.0]
        simplices[:, 3, 4] = [4.0, 0.0, 1.0]

        bvh = BoundingVolumeHierarchy(simplices, 5)

        @test bvh.tree isa BVHNode
        @test bvh.tree.leftChild isa BVHTree
        @test bvh.tree.rightChild isa BVHTree
    end

    @testset "Depth and Bounding Box Verification" begin

        simplices = zeros(Float64, 3, 3, 4)
        simplices[:, 1, 1] = [0.0, 0.0, 0.0]
        simplices[:, 1, 2] = [1.0, 0.0, 0.0]
        simplices[:, 1, 3] = [0.0, 1.0, 0.0]
        simplices[:, 1, 4] = [0.0, 0.0, 1.0]
        simplices[:, 2, 1] = [2.0, 0.0, 0.0]
        simplices[:, 2, 2] = [3.0, 0.0, 0.0]
        simplices[:, 2, 3] = [2.0, 1.0, 0.0]
        simplices[:, 2, 4] = [2.0, 0.0, 1.0]
        simplices[:, 3, 1] = [4.0, 0.0, 0.0]
        simplices[:, 3, 2] = [5.0, 0.0, 0.0]
        simplices[:, 3, 3] = [4.0, 1.0, 0.0]
        simplices[:, 3, 4] = [4.0, 0.0, 1.0]

        bvh = BoundingVolumeHierarchy(simplices, 5)

        @test bvh.tree.depth == 5
        @test bvh.tree.leftChild.depth == 4 || bvh.tree.leftChild isa BVHLeaf

        @test bvh.bbox[1, 2] == 5.0
        @test bvh.bbox[2, 1] == 0.0
    end

    @testset "Split Logic Consistency" begin

        # Original pair separated in Y
        # Original pair separated in Y
        sA = zeros(Float64, 3, 2, 4)
        sA[:, 1, 1] = [0.0, 0.0, 0.0]
        sA[:, 1, 2] = [1.0, 0.0, 0.0]
        sA[:, 1, 3] = [0.0, 1.0, 0.0]
        sA[:, 1, 4] = [0.0, 0.0, 1.0]
        sA[:, 2, 1] = [0.0, 2.0, 0.0]
        sA[:, 2, 2] = [1.0, 2.0, 0.0]
        sA[:, 2, 3] = [0.0, 3.0, 0.0]
        sA[:, 2, 4] = [0.0, 2.0, 1.0]

        bvhA = BoundingVolumeHierarchy(sA, 1)

        sB = copy(sA)
        sB[3, :, :] .+= 10.0
        bvhB = BoundingVolumeHierarchy(sB, 1)

        # Test shows that a one depth split doesn't discern y axis change
        @test bvhA.tree == bvhB.tree
    end

end
