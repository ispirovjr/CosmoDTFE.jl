using TetGen
using StaticArrays
using GLMakie

using .TesselationCore

points3d = [point3(@SVector rand(3)) for _ in 1:50]


coords, tets = TesselationCore.tesselate(points3d)

tets

#TODO Create Tesselation Constructor
#TODO Plot Tet -> Plot Tesselation (run each child)

fig = Figure()
ax = Axis3(fig[1, 1], title="Tesselation", aspect=:data)


i = 0
for tet in eachrow(tets)
    
    pos = coords[:,tet]    #get 3x4 coords
    posEs = ntuple(j -> point3(pos[:,j]),4)

    tetro = TesselationCore.Tetrahedron(posEs)
    
    TesselationCore.plotTet!(ax,tetro)
    
    i+=1
end

fig
plotTet(teter)


tester  = coords[:,tets[1,:]]


coords



_computeVolume(tester)