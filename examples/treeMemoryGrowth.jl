using StaticArrays
using JLD2
using Plots

import illustris_julia as il

include("./TesselationCore.jl")
import .TesselationCore

BVH = TesselationCore.BVH
point3 = TesselationCore.point3

#basePath = "../ThesisMaster/Illustris/"; #virgo
basePath = "../DTFE/Illustris3/output"; # laptop

fields = ["SubhaloMass","SubhaloCM"];
subhalos = il.groupcat.loadSubhalos(basePath,135,fields)

positions = subhalos["SubhaloCM"]
masses = subhalos["SubhaloMass"]

gap = 1
points = positions[:,1:gap:end]
ps = [point3(points[1,i], points[2,i], points[3,i]) for i in 1:size(points,2)]

depths = 1:18

sizes = []
sizBvh = []
sizTes = []
sizTet = []

for d in depths
    bvh,tes,tets = TesselationCore.standardEstimator(ps,masses,d)
    sizeLim = (Base.summarysize(bvh)+Base.summarysize(tes)+Base.summarysize(tets))/1e6
    push!(sizes,sizeLim)
    push!(sizBvh,Base.summarysize(bvh)/1e6)
    push!(sizTes,Base.summarysize(tes)/1e6)
    push!(sizTet,Base.summarysize(tets)/1e6)
    
end

plot(depths,sizes,
    xticks=depths,
    title = "Data Growth",
    xlabel="Tree Depth",
    ylabel="Size (Mb)",
    label="My Estimator",
)
plot!(depths,sizBvh,
    label="BVH Controbution",
)

plot!(depths,sizTes,
    label="Points Controbution",
)

plot!(depths,sizTet,
    label="Tetrahedron Controbution",
)


savefig("./Images/MyDataGrowth.png")

using PhaseSpaceDTFE 
depthsJob = 1:6
N = 32
width = 75000
simBox = SimBox(width,N)

sizeJob = []

for d in depthsJob
    estimator = DTFE_periodic(points', masses, d, simBox, pad=0)
    siz = Base.summarysize(estimator)/1e6
    println(siz)
    push!(sizeJob,siz)
end

plot(depthsJob .* 3,sizeJob,
    title = "Data Growth",
    xlabel="Tree Depth",
    ylabel="Size (Mb)",
    label="Job Estimator",
)

sizeJob

plot!(depthsJob .* 3,sizeJob,
    title = "Data Growth",
    xlabel="Tree Depth",
    ylabel="Size (Mb)",
    label="Job Estimator",
)

savefig("./Images/DataGrowthComp.png")