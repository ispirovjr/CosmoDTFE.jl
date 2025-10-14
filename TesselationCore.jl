module TesselationCore

export point2, point3, Tetrahedron, Triangulation3D, tesselate, x, y, z,plotTet,plotTet!

using StaticArrays
using TetGen

include("elements.jl")
using .Elements: Tetrahedron, point3

include("tesselate.jl")
using .Tesselate: tesselate

include("plotting.jl")
using .Plotting: plotTet, plotTet!

end 

