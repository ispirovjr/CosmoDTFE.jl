using TetGen
using StaticArrays
using JLD
using BenchmarkTools
using LinearAlgebra
using Plots

include("./TesselationCore.jl")
import .TesselationCore
using LsqFit, Statistics
import illustris_julia as il

BVH = TesselationCore.BVH
point3 = TesselationCore.point3


basePath = "../ThesisMaster/Illustris/";

fields = ["Masses","Coordinates","ParticleIDs"];

load135 = il.snapshot.loadSubset(basePath,135,"gas",fields)
gap = 10000
points = load135["Coordinates"][:,1:gap:end]

ps = [point3(points[1,i], points[2,i], points[3,i]) for i in 1:size(points,2)]

bvh,tes,tets = TesselationCore.standardEstimator(ps,12)




midPoint = mean(bvh.bbox,dims=2)[:,1]
