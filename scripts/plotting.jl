module Plotting

using GLMakie
using ..Elements

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





end