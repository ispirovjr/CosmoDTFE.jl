using Test
using CosmoDTFE

@testset "Tessellate" begin
    irregularPoints = [
        Point3(0.00, 0.00, 0.00),
        Point3(1.00, 0.00, 0.00),
        Point3(0.00, 1.00, 0.00),
        Point3(0.00, 0.00, 1.00),
        Point3(0.83, 0.31, 0.47),
        Point3(0.26, 0.72, 0.58),
        Point3(0.54, 0.46, 0.91),
        Point3(0.77, 0.88, 0.22),
    ]

    @testset "tessellate preserves input points and returns valid tetrahedra" begin
        points = irregularPoints
        coords, tets = tessellate(points)

        @test size(coords, 1) == 3
        @test size(tets, 2) == 4
        @test size(tets, 1) > 0
        @test all(tets .>= 1)
        @test all(tets .<= size(coords, 2))
        @test Set(points) == Set(Point3(coords[:, pointId]) for pointId in axes(coords, 2))
    end

    @testset "matrix and Point3 inputs produce the same geometry" begin
        points = irregularPoints
        pointMatrix = reduce(hcat, points)

        vectorCoords, vectorTets = tessellate(points)
        matrixCoords, matrixTets = tessellate(pointMatrix)

        @test matrixCoords ≈ vectorCoords
        @test matrixTets == vectorTets
    end

    @testset "tetrahedra have positive volume" begin
        points = irregularPoints
        coords, tets = tessellate(points)
        tessPoints = [Point3(coords[:, pointId]) for pointId in axes(coords, 2)]
        volumes = [CosmoDTFE.computeVolume(tessPoints[tet]) for tet in eachrow(tets)]

        @test all(volume -> volume > 0.0, volumes)
    end

    @testset "cube corner tessellation fills the unit cube volume" begin
        cubePoints = [
            Point3(0.0, 0.0, 0.0), Point3(1.0, 0.0, 0.0),
            Point3(0.0, 1.0, 0.0), Point3(1.0, 1.0, 0.0),
            Point3(0.0, 0.0, 1.0), Point3(1.0, 0.0, 1.0),
            Point3(0.0, 1.0, 1.0), Point3(1.0, 1.0, 1.0),
        ]
        coords, tets = tessellate(cubePoints)
        tessPoints = [Point3(coords[:, pointId]) for pointId in axes(coords, 2)]
        totalVolume = sum(CosmoDTFE.computeVolume(tessPoints[tet]) for tet in eachrow(tets))

        @test all(tets .>= 1)
        @test all(tets .<= size(coords, 2))
        @test totalVolume ≈ 1.0
    end
end
