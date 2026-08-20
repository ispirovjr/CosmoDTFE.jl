using Test
using CosmoDTFE

@testset "BVH" begin
    simplices = zeros(Float64, 3, 4, 3)

    simplices[:, 1, 1] = [0.0, 0.0, 0.0]
    simplices[:, 2, 1] = [1.0, 0.0, 0.0]
    simplices[:, 3, 1] = [0.0, 1.0, 0.0]
    simplices[:, 4, 1] = [0.0, 0.0, 1.0]

    simplices[:, 1, 2] = [2.0, 0.0, 0.0]
    simplices[:, 2, 2] = [3.0, 0.0, 0.0]
    simplices[:, 3, 2] = [2.0, 1.0, 0.0]
    simplices[:, 4, 2] = [2.0, 0.0, 1.0]

    simplices[:, 1, 3] = [0.25, 0.25, 0.25]
    simplices[:, 2, 3] = [0.75, 0.25, 0.25]
    simplices[:, 3, 3] = [0.25, 0.75, 0.25]
    simplices[:, 4, 3] = [0.25, 0.25, 0.75]

    @testset "global bounding box matches simplex extents" begin
        bvh = BoundingVolumeHierarchy(simplices, 3)

        @test bvh.bbox isa CosmoDTFE.SearchBox3D
        @test bvh.bbox[:, 1] ≈ [0.0, 0.0, 0.0]
        @test bvh.bbox[:, 2] ≈ [3.0, 1.0, 1.0]
    end

    @testset "depth zero stores all simplex ids in one leaf" begin
        bvh = BoundingVolumeHierarchy(simplices, 0)

        @test bvh.tree isa BVHLeaf
        @test bvh.tree.data == [1, 2, 3]
    end

    @testset "direct tetrahedron construction matches simplex construction" begin
        points = reshape(simplices, 3, :)
        pointVectors = [Point3(point) for point in eachcol(points)]
        tetrahedra = [
            6 10 4
            5 12 2
            8  9 1
            7 11 3
        ]
        tensorData = points[:, tetrahedra]
        tensorBoxes = CosmoDTFE.getBoxes(tensorData)
        vectorBoxes = CosmoDTFE.getBoxes([
            tensorData[:, :, tetrahedronId]
            for tetrahedronId in axes(tensorData, 3)
        ])
        directBoxes = CosmoDTFE.getBoxes(pointVectors, tetrahedra)

        @test size(tensorBoxes) == (3, 2, 3)
        @test tensorBoxes[:, :, 1] == [2.0 3.0; 0.0 1.0; 0.0 1.0]
        @test tensorBoxes[:, :, 2] == [0.25 0.75; 0.25 0.75; 0.25 0.75]
        @test tensorBoxes[:, :, 3] == [0.0 1.0; 0.0 1.0; 0.0 1.0]
        @test directBoxes == tensorBoxes
        @test vectorBoxes == tensorBoxes

        tensorBvh = BoundingVolumeHierarchy(tensorData, 3)
        directBvh = BoundingVolumeHierarchy(pointVectors, tetrahedra, 3)
        @test directBvh == tensorBvh
        queryPoints = (
            Point3(-0.5, -0.5, -0.5),
            Point3(0.25, 0.25, 0.25),
            Point3(1.5, 0.5, 0.5),
            Point3(2.5, 0.25, 0.25),
            Point3(3.5, 1.5, 1.5),
        )
        for point in queryPoints
            expectedCandidates = recursiveSearch(point, tensorBvh)
            @test recursiveSearch(point, directBvh) == expectedCandidates
        end

        box = [-1.0 4.0; -2.0 2.0; -3.0 3.0]
        @test BoundingVolumeHierarchy(pointVectors, tetrahedra, 2; box=box) ==
              BoundingVolumeHierarchy(tensorData, 2; box=box)

    end

    @testset "split-spanning boxes stay reachable from both sides" begin
        simplices = zeros(Float64, 3, 4, 1)
        simplices[:, 1, 1] = [0.25, 0.25, 0.25]
        simplices[:, 2, 1] = [0.75, 0.25, 0.25]
        simplices[:, 3, 1] = [0.25, 0.75, 0.25]
        simplices[:, 4, 1] = [0.25, 0.25, 0.75]

        box = [0.0 1.0; 0.0 1.0; 0.0 1.0]
        bvh = BoundingVolumeHierarchy(simplices, 1; box=box)

        @test 1 in recursiveSearch(Point3(0.25, 0.25, 0.25), bvh)
        @test 1 in recursiveSearch(Point3(0.25, 0.75, 0.25), bvh)
    end

    @testset "leaf equality ignores candidate order" begin
        @test BVHLeaf([1, 2, 3]) == BVHLeaf([3, 2, 1])
        @test BVHLeaf([1, 2, 3]) != BVHLeaf([1, 2])
    end
end
