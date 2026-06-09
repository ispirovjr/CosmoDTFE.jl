using Test
using StaticArrays
using CosmoDTFE

@testset "Tessellate" begin

    @testset "tessellate with Point3 vector" begin

        points = [Point3(rand(), rand(), rand()) for _ in 1:20]

        coords, tets = tessellate(points)

        @test size(coords, 1) == 3  # 3D coordinates
        @test size(coords, 2) >= 20  # At least as many points as input
        @test size(tets, 2) == 4  # Each tetrahedron has 4 vertices
        @test size(tets, 1) > 0  # At least one tetrahedron

        # All indices should be valid
        @test all(tets .>= 1)
        @test all(tets .<= size(coords, 2))
    end

    @testset "tessellate with Matrix" begin
        # Create random point cloud as matrix [3 x N]
        points = rand(3, 20)

        coords, tets = tessellate(points)

        @test size(coords, 1) == 3
        @test size(tets, 2) == 4
        @test size(tets, 1) > 0
    end

    @testset "tessellate preserves convex hull" begin
        # Create cube corners - should produce valid tessellation
        cubePoints = [
            Point3(0.0, 0.0, 0.0), Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0), Point3(1.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0), Point3(1.0, 0.0, 1.0),
            Point3(0.0, 1.0, 1.0), Point3(1.0, 1.0, 1.0)
        ]

        coords, tets = tessellate(cubePoints)

        @test size(tets, 1) >= 5
    end

end
