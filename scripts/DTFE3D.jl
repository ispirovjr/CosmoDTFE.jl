using StaticArrays
using GLMakie
using JLD2

using Statistics
using ColorSchemes

@load "./saves/3DdensEight.jld2" dens

normdata = dens ./ median(dens);

lowColor = get(ColorSchemes.acton, LinRange(0,1,256))[1]


fig = Figure(backgroundcolor=lowColor)
ax = LScene(fig[1,1], scenekw=(show_axis=false, backgroundcolor=lowColor))

maximum(normdata)

volume!(
    ax,
    normdata;
    algorithm  = :mip,
    colormap   = :acton,
)

fig


using JLD
JLD.save("./Images/ImprovedParticles.png", fig)