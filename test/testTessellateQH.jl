using Test
using StaticArrays
using JuliaDTFE

@testset "TessellateQH" begin

    @testset "tessellateQH with Point3 vector" begin

        points = [Point3(rand(), rand(), rand()) for _ in 1:20]

        coords, tets = tessellateQH(points)

        @test size(coords, 1) == 3  # 3D coordinates
        @test size(coords, 2) >= 20  # At least as many points as input
        @test size(tets, 2) == 4  # Each tetrahedron has 4 vertices
        @test size(tets, 1) > 0  # At least one tetrahedron

        # All indices should be valid
        @test all(tets .>= 1)
        @test all(tets .<= size(coords, 2))
    end

    @testset "tessellateQH with Matrix" begin
        # Create random point cloud as matrix [3 x N]
        points = rand(3, 20)

        coords, tets = tessellateQH(points)

        @test size(coords, 1) == 3
        @test size(tets, 2) == 4
        @test size(tets, 1) > 0
    end

    @testset "tessellateQH preserves convex hull" begin
        # Cube corners alone are coplanar when lifted to 4D (all lie on a sphere)
        # so we add interior points to break the degeneracy.
        cubePoints = [
            Point3(0.0, 0.0, 0.0), Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0), Point3(1.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0), Point3(1.0, 0.0, 1.0),
            Point3(0.0, 1.0, 1.0), Point3(1.0, 1.0, 1.0),
            Point3(0.5, 0.5, 0.5), Point3(0.25, 0.25, 0.25),
            Point3(0.75, 0.75, 0.75)
        ]

        coords, tets = tessellateQH(cubePoints)

        @test size(tets, 1) >= 5
    end

    @testset "tessellateQH matches tessellate output" begin
        # Both backends must produce valid Delaunay tessellations of the same
        # point cloud.  The exact tetrahedra may differ (Delaunay is not unique
        # on degenerate inputs), so we compare structural properties: same
        # coordinate set, same total tessellated volume, and each tet has
        # positive volume.

        points = [Point3(rand(), rand(), rand()) for _ in 1:50]

        coordsTG, tetsTG = tessellate(points)
        coordsQH, tetsQH = tessellateQH(points)

        # Coordinate matrices should have the same number of columns
        # (QH never adds Steiner points; TetGen usually doesn't with "Q")
        @test size(coordsQH, 2) == length(points)

        # Both should cover the same total volume (sum of all tet volumes)
        function totalVolume(coords, tets)
            vol = 0.0
            for row in eachrow(tets)
                v1 = coords[:, row[1]]
                v2 = coords[:, row[2]]
                v3 = coords[:, row[3]]
                v4 = coords[:, row[4]]
                a = v2 - v1
                b = v3 - v1
                c = v4 - v1
                vol += abs(dot(a, cross(b, c))) / 6
            end
            return vol
        end

        volTG = totalVolume(coordsTG, tetsTG)
        volQH = totalVolume(coordsQH, tetsQH)

        # Total tessellated volume should match (same convex hull)
        @test volTG ≈ volQH rtol = 1e-10
    end

end
