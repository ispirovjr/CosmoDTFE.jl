using Test
using CosmoDTFE
using Statistics

@testset "Estimators" begin
    function weightedMassFromVertexDensities(points, tetrahedra, rhoStar)
        return sum(eachrow(tetrahedra)) do tet
            vol = CosmoDTFE.computeVolume(points[tet])
            mean(rhoStar[tet]) * vol
        end
    end

    function sphericalCoordinateCloud()
        rMin, rMax = 0.25, 1.00
        nR, nTheta, nPhi = 20, 24, 28
        thetaMin, thetaMax = 0.15 * pi, 0.85 * pi
        dPhi = 2pi / nPhi

        points = Point3[]
        for r in range(rMin, rMax, length=nR)
            for theta in range(thetaMin, thetaMax, length=nTheta)
                for phi in range(0.0, 2pi - dPhi, length=nPhi)
                    push!(
                        points,
                        Point3(
                            r * sin(theta) * cos(phi),
                            r * sin(theta) * sin(phi),
                            r * cos(theta),
                        ),
                    )
                end
            end
        end

        return points
    end

    unitTetPoints = [
        Point3(0.0, 0.0, 0.0),
        Point3(1.0, 0.0, 0.0),
        Point3(0.0, 1.0, 0.0),
        Point3(0.0, 0.0, 1.0),
    ]
    unitTet = [1 2 3 4]

    irregularPoints = [
        Point3(0.00, 0.00, 0.00),
        Point3(1.00, 0.00, 0.00),
        Point3(0.00, 1.00, 0.00),
        Point3(0.00, 0.00, 1.00),
        Point3(1.00, 1.00, 0.00),
        Point3(1.00, 0.00, 1.00),
        Point3(0.00, 1.00, 1.00),
        Point3(1.00, 1.00, 1.00),
        Point3(0.32, 0.41, 0.27),
        Point3(0.71, 0.28, 0.54),
        Point3(0.24, 0.76, 0.61),
        Point3(0.58, 0.63, 0.84),
    ]
    irregularWeights = collect(range(0.5, 1.6, length=length(irregularPoints)))

    @testset "Triangulation density normalization is exact on one tetrahedron" begin
        weights = [1.0, 2.0, 4.0, 8.0]

        triangulation = Triangulation3D(unitTetPoints, unitTet, weights)
        volume = CosmoDTFE.computeVolume(unitTetPoints)

        @test triangulation.rhoStar ≈ 4.0 .* weights ./ volume
        @test weightedMassFromVertexDensities(unitTetPoints, unitTet, triangulation.rhoStar) ≈ sum(weights)
    end

    @testset "Density normalization is translation and scale invariant" begin
        weights = [1.0, 2.0, 4.0, 8.0]
        baseTriangulation = Triangulation3D(unitTetPoints, unitTet, weights)

        offset = Point3(-3.0, 2.0, 5.0)
        scale = 2.5
        transformedPoints = [offset + scale * point for point in unitTetPoints]
        transformedTriangulation = Triangulation3D(transformedPoints, unitTet, weights)

        @test transformedTriangulation.rhoStar ≈ baseTriangulation.rhoStar ./ scale^3
        @test weightedMassFromVertexDensities(transformedPoints, unitTet, transformedTriangulation.rhoStar) ≈ sum(weights)
    end

    @testset "dtfe performs exact barycentric interpolation inside one tetrahedron" begin
        weights = [1.0, 2.0, 3.0, 4.0]
        triangulation = Triangulation3D(unitTetPoints, unitTet, weights)
        bvh = BoundingVolumeHierarchy(reduce(hcat, unitTetPoints)[:, unitTet], 2)

        vertexDensities = triangulation.rhoStar
        lambda = [0.1, 0.2, 0.3, 0.4]
        queryPoint = Point3(0.2, 0.3, 0.4)
        facePoint = Point3(1 / 3, 1 / 3, 1 / 3)

        @test dtfe(unitTetPoints[1], bvh, unitTet, triangulation) ≈ vertexDensities[1]
        @test dtfe(queryPoint, bvh, unitTet, triangulation) ≈ sum(lambda .* vertexDensities)
        @test dtfe(facePoint, bvh, unitTet, triangulation) ≈ mean(vertexDensities[2:4])
        @test dtfe(Point3(1.01, 0.01, 0.01), bvh, unitTet, triangulation) == 0.0
    end

    @testset "DensityEstimator respects affine density scaling" begin
        estimator = DensityEstimator(irregularPoints, irregularWeights; depth=5)

        offset = Point3(4.0, -2.0, 3.0)
        scale = 3.0
        scaledPoints = [offset + scale * point for point in irregularPoints]
        scaledEstimator = DensityEstimator(scaledPoints, irregularWeights; depth=5)

        queryPoints = [
            Point3(0.25, 0.25, 0.25),
            Point3(0.45, 0.37, 0.49),
            Point3(0.72, 0.61, 0.58),
        ]

        for queryPoint in queryPoints
            @test scaledEstimator(offset + scale * queryPoint) ≈ estimator(queryPoint) / scale^3 
        end
    end

    @testset "DensityEstimator conserves total weighted mass" begin
        estimator = DensityEstimator(irregularPoints, irregularWeights; depth=5)

        totalMass = weightedMassFromVertexDensities(estimator.triangulation.points, estimator.tetrahedra, estimator.triangulation.rhoStar)

        @test totalMass ≈ sum(irregularWeights)
    end

    @testset "DensityEstimator recovers spherical coordinate sampling profile" begin
        points = sphericalCoordinateCloud()
        estimator = DensityEstimator(points; depth=8)

        radii = collect(range(0.40, 0.92, length=7))
        thetaSamples = collect(range(0.30pi, 0.70pi, length=9))
        phiSamples = collect(range(0.0, 2pi - 2pi / 16, length=16))
        radialProfile = [
            mean(
                estimator(Point3(
                    radius * sin(theta) * cos(phi),
                    radius * sin(theta) * sin(phi),
                    radius * cos(theta),
                )) * sin(theta)
                for theta in thetaSamples for phi in phiSamples
            )
            for radius in radii
        ]
        expectedRatios = [(radii[i + 1] / radii[i])^2 for i in 1:length(radii)-1]
        measuredRatios = [radialProfile[i] / radialProfile[i + 1] for i in 1:length(radialProfile)-1]

        @test all(isfinite, radialProfile)
        @test all(density -> density > 0.0, radialProfile)
        @test issorted(radialProfile; rev=true)

        for ratioId in eachindex(expectedRatios)
            @test measuredRatios[ratioId] ≈ expectedRatios[ratioId] rtol = 0.08
        end
    end

    @testset "DensityEstimator returns zero outside the tessellated domain" begin
        estimator = DensityEstimator(irregularPoints; depth=5)

        @test estimator(Point3(100.0, 100.0, 100.0)) == 0.0
        @test estimator(Point3(-1.0, -1.0, -1.0)) == 0.0
    end

    
end
