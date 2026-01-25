using Test
using JuliaDTFE

@testset "JuliaDTFE" begin
    include("testElements.jl")
    include("testBvh.jl")
    include("testTessellate.jl")
    include("testSearchers.jl")
    include("testEstimators.jl")
end

