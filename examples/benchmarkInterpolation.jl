using StaticArrays
using JLD2
using LinearAlgebra
using Plots
using BenchmarkTools


include("./TesselationCore.jl")
import .TesselationCore

BVH = TesselationCore.BVH
point3 = TesselationCore.point3


@load "./saves/dtfeEstimatorTenth.jld2" bvh tes tets


idRand = 521
tet = tets[idRand,:]
simp = tes.points[tet]
rhos = tes.ρStar[tet]



@inline function invertClassic(rhos,simplex) 
    r = rhos[2:end] .- rhos[1]

    v1, v2, v3, v4 = simplex[1], simplex[2], simplex[3], simplex[4]

    a = v2 - v1 
    b = v3 - v1 
    c = v4 - v1

    mat = SMatrix{3,3}(hcat(a, b, c))    
    return inv(mat)*r #
end

function interpOld(rhos, simp,point)
    delRho=invertClassic(rhos,simp)
    interpolation = rhos[1] + dot((point - simp[1]),delRho)
    return interpolation
end

function interpNew(rhos,simp,point)
    v1 = simp[1]
    invM = inv(SMatrix{3,3}(hcat(simp[2]-v1, simp[3]-v1, simp[4]-v1)))
    diff = point - v1
    λ234 = invM * diff
    λ1   = 1 - sum(λ234)
    interpolation = λ1*rhos[1] + λ234[1]*rhos[2] + λ234[2]*rhos[3] + λ234[3]*rhos[4]
    return interpolation

end

simp

pt = mean(simp,dims=1)[1]


simp[1]

interpOld(rhos,simp,pt)
interpNew(rhos,simp,pt)

test1 = @benchmarkable interpOld(rhos,simp,pt)
test2 = @benchmarkable interpNew(rhos,simp,pt)

stats1 = run(test1)
stats2 = run(test2)


mins = mapreduce(x -> x, (a,b) -> min.(a,b), simp)
maxs = mapreduce(x -> x, (a,b) -> max.(a,b), simp)

N = 10

xs = range(mins[1], maxs[1], length=N)
ys = range(mins[2], maxs[2], length=N)
zs = range(mins[3], maxs[3], length=N)

valsOld = [interpOld(rhos,simp,[x,y,z]) for x in xs, y in ys, z in zs];

hmold = heatmap(xs,ys,valsOld[:,:,5],title = "Old Interpolation")
for i in 1:4
    if i ==4
        plot!([simp[4][1],simp[1][1]],[simp[4][2],simp[1][2]],label="Tetrahedron",color=:black)
        break
    end
    plot!([simp[i][1],simp[i+1][1]],[simp[i][2],simp[i+1][2]],label="",color=:black)

end

hmold

valsNew = [interpNew(rhos,simp,[x,y,z]) for x in xs, y in ys, z in zs];



hmnew = heatmap(xs,ys,sum(valsNew,dims=3)[:,:,1],title = "New Interpolation")
for i in 1:4
    if i ==4
        plot!([simp[4][1],simp[1][1]],[simp[4][2],simp[1][2]],label="Tetrahedron",color=:black)
        break
    end
    plot!([simp[i][1],simp[i+1][1]],[simp[i][2],simp[i+1][2]],label="",color=:black)

end

hmnew

plot(hmold, hmnew, layout=(1,2),size=(1000,400))
savefig("./Images/interpolationComp.png")
