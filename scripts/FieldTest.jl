using TetGen
using StaticArrays
using JLD
using BenchmarkTools
using LinearAlgebra
using Plots

import illustris_julia as il

include("./TesselationCore.jl")
import .TesselationCore

BVH = TesselationCore.BVH
point3 = TesselationCore.point3


basePath = "../ThesisMaster/Illustris/";

fields = ["SubhaloMass","SubhaloCM"];
subhalos = il.groupcat.loadSubhalos(basePath,135,fields)

positions = subhalos["SubhaloCM"]

gap = 50
points = positions[:,1:gap:end]
ps = [point3(points[1,i], points[2,i], points[3,i]) for i in 1:size(points,2)]

bvh,tes,tets = TesselationCore.standardEstimator(ps,10)

N = 1000

width = 75000

step = width/N


xs = bvh.bbox[1,1]:step:bvh.bbox[1,2]
ys = bvh.bbox[2,1]:step:bvh.bbox[2,2]

z = (bvh.bbox[3,2] + bvh.bbox[3,1])/2

dens = zeros(N,N)

println("Slice")

for (i,x) in pairs(xs)
    for (j,y) in pairs(ys)
        dens[i,j] = TesselationCore.DTFE([x,y,z],bvh,tets,tes)

    end
end

med = median(dens)

Plots.heatmap(dens ./med,clim=(0,25))

savefig("./Images/DenSlice.png")


xs = bvh.bbox[1,1]:step:bvh.bbox[1,2]
ys = bvh.bbox[2,1]:step:bvh.bbox[2,2]

zs =  bvh.bbox[3,1]:step:bvh.bbox[3,2]

dens = zeros(N,N,N)

println("Fat Chunk")

for (i,x) in pairs(xs)
    for (j,y) in pairs(ys)
        for (k,z) in pairs(zs)
            dens[i,j,k] = TesselationCore.DTFE([x,y,z],bvh,tets,tes)

        end

    end
end



save("../saves/3Ddens.jld", dens)  