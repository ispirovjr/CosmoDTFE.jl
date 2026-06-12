using Test
using CosmoDTFE

"""
This file runs all the tests for the CosmoDTFE package.

Elements - tests the basic data structures and functions.
Bvh - tests the bounding volume hierarchy.
Tessellate - tests the tessellation algorithm.
Searchers - tests the search algorithms.
Estimators - tests the complete density estimators.
"""


@testset "CosmoDTFE" begin
    include("testElements.jl")
    include("testBvh.jl")
    include("testTessellate.jl")
    include("testSearchers.jl")
    include("testEstimators.jl")
    include("testVelocity.jl")
    include("testPeriodic.jl")
    include("testPhaseSpace.jl")
    include("testComposite.jl")
    if get(ENV, "DO_TETGEN_THREAD_TESTS", "false") == "true"
        include("tetgenThreadSafety.jl")
    else
        @info "Skipping TetGen thread safety test. Set DO_TETGEN_THREAD_TESTS=true to run."
    end
end

