using Test
using StaticArrays
using CosmoDTFE

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
    separatedTets = [1 2 3 4; 5 6 7 8]

    @testset "intersection3D handles interiors and boundaries" begin
        simplex = SMatrix{4,3,Float64}([
            0.0 0.0 0.0;
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])

        reversedSimplex = SMatrix{4,3,Float64}([
            0.0 0.0 1.0;
            0.0 1.0 0.0;
            1.0 0.0 0.0;
            0.0 0.0 0.0
        ])

        @test CosmoDTFE.intersection3D(Point3(0.25, 0.25, 0.25), simplex)
        @test CosmoDTFE.intersection3D(Point3(0.25, 0.25, 0.25), reversedSimplex)
        @test CosmoDTFE.intersection3D(Point3(0.0, 0.0, 0.0), simplex)
        @test CosmoDTFE.intersection3D(Point3(1/3, 1/3, 1/3), simplex)
        @test !CosmoDTFE.intersection3D(Point3(0.34, 0.34, 0.34), simplex)
        @test !CosmoDTFE.intersection3D(Point3(-0.01, 0.25, 0.25), simplex)
    end

    @testset "findId returns deterministic simplex ids" begin
        simplices = separatedTetPoints[separatedTets]
        bvh = BoundingVolumeHierarchy(reduce(hcat, separatedTetPoints)[:, separatedTets], 4)

        @test findId(Point3(0.25, 0.25, 0.25), simplices, bvh) == 1
        @test findId(Point3(2.25, 0.25, 0.25), simplices, bvh) == 2
        @test findId(Point3(10.0, 10.0, 10.0), simplices, bvh) === nothing
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
        tetrahedra = [1 2 3 4; 1 2 3 4]
        bvh = BoundingVolumeHierarchy(reduce(hcat, points)[:, tetrahedra], 4)

        @test findAllIds(Point3(0.25, 0.25, 0.25), points, tetrahedra, bvh) == [1, 2]
        @test isempty(findAllIds(Point3(2.0, 2.0, 2.0), points, tetrahedra, bvh))
    end
end
