using Test
using StaticArrays
using JuliaDTFE

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

    @testset "Tetrahedron construction" begin
        # Create a simple tetrahedron
        v1 = Point3(0.0, 0.0, 0.0)
        v2 = Point3(1.0, 0.0, 0.0)
        v3 = Point3(0.0, 1.0, 0.0)
        v4 = Point3(0.0, 0.0, 1.0)
        verts = (v1, v2, v3, v4)

        tet = Tetrahedron(verts)

        @test tet.verts == verts
        @test tet.vol ≈ 1 / 6  # Volume of unit tetrahedron
    end

    @testset "Triangulation3D construction" begin
        # TODO: Add tests for Triangulation3D with tessellation
        # This requires integration with TetGen
        @test_skip "Triangulation3D requires tessellation integration"
    end

end
