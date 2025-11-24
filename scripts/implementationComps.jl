using StaticArrays
using JLD2
using BenchmarkTools
using LinearAlgebra
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

bvh,tes,tets = TesselationCore.standardEstimator(ps,masses,9) 

@save "./saves/dtfeEstimator.jld2" bvh tes tets # 25.7 Mb for depth 9

print((Base.summarysize(bvh)+Base.summarysize(tes)+Base.summarysize(tets))/1e6)


N = 32
width = 75000
step = width/N

xs = bvh.bbox[1,1]:step:bvh.bbox[1,2]
ys = bvh.bbox[2,1]:step:bvh.bbox[2,2]
zs =  bvh.bbox[3,1]:step:bvh.bbox[3,2]

Threads.nthreads()

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 900
BenchmarkTools.DEFAULT_PARAMETERS.samples = 5

linBench = @benchmarkable [TesselationCore.DTFE([x, y, z], bvh, tets, tes) for x in xs, y in ys,z in zs]
parallelBench = @benchmarkable TesselationCore.DTFEMultiThread([xs,ys,zs],bvh,tets,tes)

parallelStat = run(parallelBench)
linStat = run(linBench)

mean(linStat)
std(linStat)
#235.360±1.345s
#16.17 GiB

#124.525±0.5 s --> single core on norma
#16.17 GiB

mean(parallelStat)
std(parallelStat)
#67.923±9.346s
#16.17 GiB

#6.175±0.1s --> 40 cores on norma
#16.17 GiB



@load "./saves/dtfeEstimator.jld2" bvh tes tets



using PhaseSpaceDTFE # run on norma - this is impossible locally

simBox = SimBox(width,N)
estimator = DTFE_periodic(points', masses, 3, simBox, pad=0)
psBench = @benchmarkable [PhaseSpaceDTFE.density([x, y, z], estimator) for x in xs, y in ys, z in zs]

@save "./saves/psEst.jld2" estimator # 121 Mb


psStat = run(psBench)

mean(psStat)
std(psStat)




