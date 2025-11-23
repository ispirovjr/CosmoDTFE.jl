module Plotting

using GLMakie
using ..Elements
using JLD
using ColorSchemes
using Statistics
using GLFW

export plotTet


function plotTet(tet::Tetrahedron; color=:blue, linewidth=2)
    fig = Figure()
    ax = Axis3(fig[1, 1], title="Tetrahedron", aspect=:data)

    v = tet.verts

    edges = [
        (1, 2), (1, 3), (1, 4),
        (2, 3), (2, 4), (3, 4)
    ]

    for (i, j) in edges
        x1, y1, z1 = v[i]
        x2, y2, z2 = v[j]
        lines!(ax, [x1, x2], [y1, y2], [z1, z2];
               color=color, linewidth=linewidth)
    end

    xs,ys,zs = [vi[1] for vi in v],[vi[2] for vi in v],[vi[3] for vi in v]
    GLMakie.scatter!(ax, xs,ys,zs; color=:red, markersize=10)

    fig
    return fig
end

function plotTet(tet::Tetrahedron; color=:blue, linewidth=2)
    fig = Figure()
    ax = Axis3(fig[1, 1], title="Tetrahedron", aspect=:data)

    v = tet.verts

    edges = [
        (1, 2), (1, 3), (1, 4),
        (2, 3), (2, 4), (3, 4)
    ]

    for (i, j) in edges
        x1, y1, z1 = v[i]
        x2, y2, z2 = v[j]
        lines!(ax, [x1, x2], [y1, y2], [z1, z2];
               color=color, linewidth=linewidth)
    end

    xs,ys,zs = [vi[1] for vi in v],[vi[2] for vi in v],[vi[3] for vi in v]
    GLMakie.scatter!(ax, xs,ys,zs; color=:red, markersize=10)

    fig
    return fig
end

function plotTet!(ax::Axis3,tet::Tetrahedron; color=:blue, linewidth=2)
     v = tet.verts

    edges = [
        (1, 2), (1, 3), (1, 4),
        (2, 3), (2, 4), (3, 4)
    ]

    for (i, j) in edges
        x1, y1, z1 = v[i]
        x2, y2, z2 = v[j]
        lines!(ax, [x1, x2], [y1, y2], [z1, z2];
               color=color, linewidth=linewidth)
    end

    xs,ys,zs = [vi[1] for vi in v],[vi[2] for vi in v],[vi[3] for vi in v]
    GLMakie.scatter!(ax, xs,ys,zs; color=:red, markersize=10)

end

"""
    render_volume(jldfile; dataset="dens", outfile="volume.png")

Load a 3D array from a JLD/JLD2 file and render it as a volume plot to `outfile`.
Works headlessly (e.g. SSH without screen).
"""
function renderVolume(jldfile; dataset="dens", outfile="volume.png")

    # --- Load array ---
    @assert isfile(jldfile) "JLD file not found: $jldfile"
    @info "Loading $dataset from $jldfile"
    data = load(jldfile, dataset)

    @assert ndims(data) == 3 "Dataset must be 3D."

    # --- Normalize data ---
    normdata = data ./ median(data)

    # --- Background color from Acton colormap ---
    lowColor = get(ColorSchemes.acton, LinRange(0,1,256))[1]


    GLMakie.activate!()    # ensure GL backend
    GLFW.WindowHint(GLFW.VISIBLE, false)   # no on-screen window

    fig = Figure(backgroundcolor=lowColor)
    ax = LScene(fig[1,1], scenekw=(show_axis=false, backgroundcolor=lowColor))

    # Volume plot
    volume!(
        ax,
        normdata;
        algorithm  = :mip,
        colormap   = :acton,
        colorrange = (0.0, 25)
    )


    # --- Save output ---
    save(outfile, fig)
    @info "Saved volume render to $outfile"
end

renderVolume("./saves/3Ddens.jld"; dataset="Density", outfile="./Images/densRender.png")

#data = load("./saves/3Ddens.jld")

end