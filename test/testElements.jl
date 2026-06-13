using Test
using StaticArrays
using CosmoDTFE

@testset "Elements" begin

    @testset "Point3 type" begin
        # Test Point3 is correctly aliased to SVector{3, Float64}
        p = Point3(1.0, 2.0, 3.0)
        @test p isa SVector{3,Float64}
        @test p[1] == 1.0
        @test p[2] == 2.0
        @test p[3] == 3.0

        # Test arithmetic operations
        p2 = Point3(4.0, 5.0, 6.0)
        @test p + p2 == Point3(5.0, 7.0, 9.0)
    end

    @testset "computeVolume correctness" begin

        v1 = Point3(0.0, 0.0, 0.0)
        v2 = Point3(1.0, 0.0, 0.0)
        v3 = Point3(0.0, 1.0, 0.0)
        v4 = Point3(0.0, 0.0, 1.0)
        verts = (v1, v2, v3, v4)

        @test CosmoDTFE.computeVolume(verts) ≈ 1 / 6

        # Scaled tetrahedron: scale by 2 -> volume scales by 8
        scaledVerts = (2 * v1, 2 * v2, 2 * v3, 2 * v4)
        @test CosmoDTFE.computeVolume(scaledVerts) ≈ 8 / 6

        # Matrix input
        vertMatrix = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
        @test CosmoDTFE.computeVolume(vertMatrix) ≈ 1 / 6

        reversedVerts = (v4, v3, v2, v1)
        @test CosmoDTFE.computeVolume(reversedVerts) ≈ 1 / 6

        degenerateVerts = (v1, v2, v3, Point3(0.25, 0.25, 0.0))
        @test CosmoDTFE.computeVolume(degenerateVerts) ≈ 0.0
    end

    @testset "Triangulation3D construction" begin
        points = [
            Point3(0.00, 0.00, 0.00),
            Point3(1.00, 0.00, 0.00),
            Point3(0.00, 1.00, 0.00),
            Point3(0.00, 0.00, 1.00),
            Point3(0.83, 0.31, 0.47),
            Point3(0.26, 0.72, 0.58),
            Point3(0.54, 0.46, 0.91),
            Point3(0.77, 0.88, 0.22),
        ]
        coords, tets = tessellate(points)

        triangulation = Triangulation3D(points, tets)

        @test length(triangulation.points) == length(points)
        @test length(triangulation.rhoStar) == length(points)
        @test all(isfinite.(triangulation.rhoStar))
        @test all(triangulation.rhoStar .> 0)  # Densities should be positive

        # Bounds check tests
        bad_weights = ones(length(points) - 1)
        @test_throws ArgumentError Triangulation3D(points, tets, bad_weights)
    end

end
