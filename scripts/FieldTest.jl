using TetGen
using StaticArrays
using JLD
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

gap = 1
points = positions[:,1:gap:end]
ps = [point3(points[1,i], points[2,i], points[3,i]) for i in 1:size(points,2)]

bvh,tes,tets = TesselationCore.standardEstimator(ps,10)

N = 512

width = 75000

step = width/N


xs = bvh.bbox[1,1]:step:bvh.bbox[1,2]
ys = bvh.bbox[2,1]:step:bvh.bbox[2,2]

z = (bvh.bbox[3,2] + bvh.bbox[3,1])/2


println("Slice")

dens = TesselationCore.DTFEMultiThread([xs,ys,z],bvh,tets,tes)[:,:,1]
med = median(dens)

Plots.heatmap(dens ./med,clim=(0,25))

savefig("./Images/DenSlice.png")


nSmall = 32
dist = step*nSmall/2

zs = z-dist:step:z+dist

dens = TesselationCore.DTFEMultiThread([xs,ys,zs],bvh,tets,tes)
med = median(dens)

den = mean(dens ./med,dims=3)[:,:,1]

Plots.heatmap(den,clim=(0,10))

savefig("./Images/DenChunk.png")



N = 256
width = 75000
step = width/N

xs = bvh.bbox[1,1]:step:bvh.bbox[1,2]
ys = bvh.bbox[2,1]:step:bvh.bbox[2,2]
zs =  bvh.bbox[3,1]:step:bvh.bbox[3,2]


println("Fat Chunk")

dens = TesselationCore.DTFEMultiThread([xs,ys,zs],bvh,tets,tes)


save("../saves/3Ddens.jld", dens, )  