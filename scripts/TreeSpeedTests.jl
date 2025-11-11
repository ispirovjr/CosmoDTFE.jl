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


basePath = "./../DTFE/Illustris3/output";

fields = ["Masses","Coordinates","ParticleIDs"];

load135 = il.snapshot.loadSubset(basePath,135,"gas",fields)
gap = 10000
points = load135["Coordinates"][:,1:gap:end]

load135 = nothing
ps = [point3(points[1,i], points[2,i], points[3,i]) for i in 1:size(points,2)]

bvh,tes,tets = TesselationCore.standardEstimator(ps,12)

Ns = [1,5,9,12,15,16,17,18,19,20,21,22,23,24]
Ts = []
errs = []
statses = []

for n in Ns

    bench = @benchmarkable TesselationCore.standardEstimator(ps,$n)
    stats = run(bench)
    push!(Ts,mean(stats).time)
    push!(errs,std(stats).time)
    push!(statses,stats)
    
    println(n)
end

model(N, p) = p[1] .* (2 .^N) .+ p[2]

p0 = [1.0, 1e8]  # initial guess for [a, b]



fit = LsqFit.curve_fit(model, Ns, Ts, p0,inplace=false)

params = coef(fit)
errors = stderror(fit)

Plots.plot(Ns,Ts,yerr = errs,yscale=:log10,label="Data")
Plots.plot!(Ns,model(Ns,params),yscale=:log10,label = "Exponential Model")
savefig("./Images/TreeGrowth.png")

means = []
meds = []

results = []

midPoint = mean(bvh.bbox,dims=2)[:,1]

for n in Ns

    bvh,tes,tets = TesselationCore.standardEstimator(ps,n)


    bench = @benchmarkable TesselationCore.DTFE(midPoint,bvh,tets,tes)
    stats = run(bench)
    push!(meds,median(stats).time)
    push!(means,mean(stats).time)
    
    println(n)
end

Plots.plot(Ns,meds,label="Medians")
Plots.plot!(Ns,means,label="Means")
savefig("./Images/EstimationTime.png")

function countLeaves(tree)
    count = 0
    
    if typeof(tree.leftChild) == TesselationCore.Bvh.BVHLeaf
        count +=1
    else 
        count += countLeaves(tree.leftChild)
    end
    if typeof(tree.rightChild) == TesselationCore.Bvh.BVHLeaf
        count +=1
    else 
        count += countLeaves(tree.rightChild)
    end

    return count

end

leavz = []

Ns = [1,5,9,12,15,16,17,18,19,20]

for n in Ns

    bvh,tes,tets = TesselationCore.standardEstimator(ps,n)

    push!(leavz,countLeaves(bvh.tree))
    
    println(n)
end

Plots.plot(Ns,leavz - 2 .^Ns)
savefig("./Images/Leaves.png")