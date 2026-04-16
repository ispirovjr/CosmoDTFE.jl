using Test
using JuliaDTFE
using StaticArrays

@testset "TetGen Thread Safety" begin
    nSets = 8
    nPoints = 10_000

    sharedPoints = [SVector{3,Float64}(rand(), rand(), rand()) for _ in 1:nPoints]

    coordResults = Vector{Any}(undef, nSets)
    tetResults = Vector{Any}(undef, nSets)

    Threads.@threads for i in 1:nSets
        coords, tets = tessellate(copy(sharedPoints))
        coordResults[i] = coords
        tetResults[i] = tets
    end

    @testset "All tessellations complete" begin
        @test all(r -> r isa AbstractMatrix, coordResults)
        @test all(r -> r isa AbstractMatrix, tetResults)
    end

    @testset "Deterministic results across threads" begin
        @test all(c -> c == coordResults[1], coordResults)
        @test all(t -> t == tetResults[1], tetResults)
    end
end
