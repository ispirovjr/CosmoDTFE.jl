module TesselationCore

export Triangulation3D, tesselate, BVH, findSimplex, findID, point3

using StaticArrays
using TetGen

include("elements.jl")
using .Elements: point3, Triangulation3D

include("tesselate.jl")
using .Tesselate: tesselate

include("bvh.jl")
using .Bvh: BVH, BVHTree, BVHLeaf, BVHNode

include("searchers.jl")
using .Searchers: recursiveSearch,  findSimplex, findID


end 

