using Test
using StaticArrays
using CosmoDTFE

searchAllocations(point, bvh) = @allocated recursiveSearch(point, bvh)

function referenceRecursiveSearch(point, leaf::BVHLeaf, searchBox::Matrix)
    return leaf.data
end

function referenceRecursiveSearch(point, tree::BVHNode, searchBox::Matrix)
    ax = tree.depth % 3 + 1
    newBox = copy(searchBox)
    line = (searchBox[ax, 2] + searchBox[ax, 1]) / 2

    if point[ax] < line
        newBox[ax, 2] = line
        return referenceRecursiveSearch(point, tree.leftChild, newBox)
    end

    newBox[ax, 1] = line
    return referenceRecursiveSearch(point, tree.rightChild, newBox)
end

@testset "Searchers" begin
    separatedTetPoints = [
        Point3(0.0, 0.0, 0.0),
        Point3(1.0, 0.0, 0.0),
        Point3(0.0, 1.0, 0.0),
        Point3(0.0, 0.0, 1.0),
        Point3(2.0, 0.0, 0.0),
        Point3(3.0, 0.0, 0.0),
        Point3(2.0, 1.0, 0.0),
        Point3(2.0, 0.0, 1.0),
    ]
    separatedTets = [1 5; 2 6; 3 7; 4 8]

    @testset "BVH traversal is allocation free and preserves candidate routing" begin
        traversalPoints = vcat(
            separatedTetPoints,
            [
                Point3(0.0, 0.0, 0.0),
                Point3(3.0, 0.0, 0.0),
                Point3(0.0, 1.0, 0.0),
                Point3(0.0, 0.0, 1.0),
            ],
        )
        traversalTets = hcat(separatedTets, [9, 10, 11, 12])
        traversalData = reduce(hcat, traversalPoints)[:, traversalTets]
        bvh = BoundingVolumeHierarchy(traversalData, 6)
        point = Point3(0.25, 0.25, 0.25)

        recursiveSearch(point, bvh)
        @test searchAllocations(point, bvh) == 0

        referenceBox = Matrix(bvh.bbox)
        queryPoints = [
            Point3(-10.0, 0.25, 0.25),
            Point3(10.0, 0.25, 0.25),
            Point3(1.5, 0.25, 0.25),
            Point3(0.15, 0.85, 0.35),
            Point3(0.95, -0.25, 1.25),
            Point3(2.05, 0.65, -0.15),
            Point3(2.85, 1.25, 0.75),
        ]

        for queryPoint in queryPoints
            expected = referenceRecursiveSearch(queryPoint, bvh.tree, referenceBox)
            @test recursiveSearch(queryPoint, bvh) === expected
        end

        @test recursiveSearch(Point3(1.5, 0.25, 0.25), bvh) == [2, 3]

        depthZero = BoundingVolumeHierarchy(
            reduce(hcat, separatedTetPoints)[:, separatedTets],
            0,
        )
        @test recursiveSearch(point, depthZero) === depthZero.tree.data
    end

    @testset "intersection3D handles interiors and boundaries" begin
        simplex = SMatrix{3,4,Float64}([
            0.0 1.0 0.0 0.0;
            0.0 0.0 1.0 0.0;
            0.0 0.0 0.0 1.0
        ])

        reversedSimplex = SMatrix{3,4,Float64}([
            0.0 0.0 1.0 0.0;
            0.0 1.0 0.0 0.0;
            1.0 0.0 0.0 0.0
        ])

        @test CosmoDTFE.intersection3D(Point3(0.25, 0.25, 0.25), simplex)
        @test CosmoDTFE.intersection3D(Point3(0.25, 0.25, 0.25), reversedSimplex)
        @test CosmoDTFE.intersection3D(Point3(0.0, 0.0, 0.0), simplex)
        @test CosmoDTFE.intersection3D(Point3(1/3, 1/3, 1/3), simplex)
        @test !CosmoDTFE.intersection3D(Point3(0.34, 0.34, 0.34), simplex)
        @test !CosmoDTFE.intersection3D(Point3(-0.01, 0.25, 0.25), simplex)
        @test CosmoDTFE.intersection3D(Point3(0.25, 0.25, 0.25), Matrix(simplex))
    end

    @testset "findId returns deterministic simplex ids" begin
        simplices = separatedTetPoints[separatedTets]
        bvh = BoundingVolumeHierarchy(reduce(hcat, separatedTetPoints)[:, separatedTets], 4)

        @test findId(Point3(0.25, 0.25, 0.25), simplices, bvh) == 1
        @test findId(Point3(2.25, 0.25, 0.25), simplices, bvh) == 2
        @test findId(Point3(10.0, 10.0, 10.0), simplices, bvh) === nothing
        @test findSimplex(Point3(0.25, 0.25, 0.25), simplices, bvh) == simplices[:, 1]
    end

    @testset "memory-efficient findId matches matrix-simplex findId" begin
        simplices = separatedTetPoints[separatedTets]
        bvh = BoundingVolumeHierarchy(reduce(hcat, separatedTetPoints)[:, separatedTets], 4)

        for queryPoint in (Point3(0.25, 0.25, 0.25), Point3(2.25, 0.25, 0.25))
            @test findId(queryPoint, separatedTetPoints, separatedTets, bvh) == findId(queryPoint, simplices, bvh)
        end
    end

    @testset "findAllIds returns every overlapping tetrahedron" begin
        points = [
            Point3(0.0, 0.0, 0.0),
            Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0),
        ]
        tetrahedra = [1 1; 2 2; 3 3; 4 4]
        bvh = BoundingVolumeHierarchy(reduce(hcat, points)[:, tetrahedra], 4)

        @test findAllIds(Point3(0.25, 0.25, 0.25), points, tetrahedra, bvh) == [1, 2]
        @test isempty(findAllIds(Point3(2.0, 2.0, 2.0), points, tetrahedra, bvh))
    end
end
