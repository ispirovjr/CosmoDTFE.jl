using Test
using StaticArrays
using LinearAlgebra
using JuliaDTFE

@testset "Searchers" begin

    @testset "intersection3D correctness" begin
        # Create a unit tetrahedron
        simplex = SMatrix{4,3,Float64}([
            0.0 0.0 0.0;
            1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
        ])

        # Point inside tetrahedron (centroid)
        centroid = SVector{3,Float64}(0.25, 0.25, 0.25)
        @test JuliaDTFE.intersection3D(centroid, simplex) == true

        # Point outside tetrahedron
        outside = SVector{3,Float64}(2.0, 2.0, 2.0)
        @test JuliaDTFE.intersection3D(outside, simplex) == false

        # Point at vertex (boundary case)
        atVertex = SVector{3,Float64}(0.0, 0.0, 0.0)
        @test JuliaDTFE.intersection3D(atVertex, simplex) == true

        # Point inside tetrahedron (centroid)

    end

    @testset "intersection3D with Vector{Float64}" begin
        # Test the Vector{Float64} overload
        simplex = [
            Point3(0.0, 0.0, 0.0),
            Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0)
        ]

        inside = [0.25, 0.25, 0.25]
        @test JuliaDTFE.intersection3D(inside, simplex) == true

        outside = [2.0, 2.0, 2.0]
        @test JuliaDTFE.intersection3D(outside, simplex) == false
    end

    @testset "findId returns correct index" begin
        # TODO: Integration test with BVH
        # This requires building a full BVH structure
        @test_skip "findId requires BVH integration test"
    end

end
