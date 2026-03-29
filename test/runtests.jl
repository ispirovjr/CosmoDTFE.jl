using Test
using JuliaDTFE

"""
This file runs all the tests for the JuliaDTFE package.

Elements - tests the basic data structures and functions.
Bvh - tests the bounding volume hierarchy.
Tessellate - tests the tessellation algorithm.
Searchers - tests the search algorithms.
Estimators - tests the complete density estimators.
"""


@testset "JuliaDTFE" begin
    include("testElements.jl")
    include("testBvh.jl")
    include("testTessellate.jl")
    include("testSearchers.jl")
    include("testEstimators.jl")
    include("testVelocity.jl")
    include("testTessellateQH.jl")
    include("testEstimatorsQH.jl")
    include("testVelocityQH.jl")
end

