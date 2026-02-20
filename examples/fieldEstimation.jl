using StaticArrays
using JLD2
using LinearAlgebra
using Plots

import illustris_julia as il

include("./TesselationCore.jl")
import .TesselationCore

BVH = TesselationCore.BVH
point3 = TesselationCore.point3


#basePath = "../ThesisMaster/Illustris/";
basePath = "../DTFE/Illustris3/output"

fields = ["SubhaloMass","SubhaloCM"];
subhalos = il.groupcat.loadSubhalos(basePath,135,fields)


positions = subhalos["SubhaloCM"]

points = positions#[:,mask]
ps = [point3(points[1,i], points[2,i], points[3,i]) for i in 1:size(points,2)]

bvh,tes,tets = TesselationCore.standardEstimator(ps,12);

@save "./saves/dtfeEstimatorTenth.jld2" bvh tes tets

N = 128

width = 75000#lim

step = width/N

xs = bvh.bbox[1,1]:step:bvh.bbox[1,2]
ys = bvh.bbox[2,1]:step:bvh.bbox[2,2]
zs =  bvh.bbox[3,1]:step:bvh.bbox[3,2]


println("Fat Chunk")

dens = TesselationCore.DTFEMultiThread([xs,ys,zs],bvh,tets,tes)

@save "./saves/3DhaloGood.jld2" dens

@load "./saves/3DdensEight.jld2" dens

using Statistics

dens = dens ./ median(dens)

minimum(dens)

Plots.heatmap(sum(dens,dims=3)[:,:,1])

Plots.heatmap(sum(dens[5:end-5,5:end-5,31:35],dims=3)[:,:,1])

