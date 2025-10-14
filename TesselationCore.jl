module TesselationCore

export point3, Tetrahedron, Triangulation3D, tesselate, x, y, z, plotTet, plotTet!,Triangulation3D

using StaticArrays
using TetGen

include("elements.jl")
using .Elements: Tetrahedron, point3, Triangulation3D

include("tesselate.jl")
using .Tesselate: tesselate

include("plotting.jl")
using .Plotting: plotTet, plotTet!

end 

