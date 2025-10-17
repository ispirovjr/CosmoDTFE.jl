using Plots
using TetGen
using StaticArrays

using .TesselationCore

points3d = [point3(@SVector rand(3)) for _ in 1:50]


coords, tets = tesselate(points3d)

Plots.scatter(coords[1,:],coords[2,:])
for tet in eachrow(tets)
    pts = hcat(coords[:,tet],coords[:,tet[1]])
    Plots.plot!(pts[1,:],pts[2,:],label="",color=:black)
end
savefig("../Images/triangulation2d.png")

