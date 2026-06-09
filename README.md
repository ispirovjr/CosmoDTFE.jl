# CosmoDTFE

[![Julia](https://img.shields.io/badge/julia-v1.9+-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

CosmoDTFE is a Julia implementation of the Delaunay Tessellation Field Estimator
(DTFE) for cosmological density and velocity fields.

DTFE estimates continuous fields from discrete samples by tessellating the point
distribution, assigning per-vertex values, and interpolating inside the
tetrahedron that contains each query point. CosmoDTFE uses TetGen for 3D
Delaunay tessellation and a bounding volume hierarchy for fast spatial lookup.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ispirovjr/CosmoDTFE.jl")
```

For local development:

```julia
] activate .
] instantiate
] test
```

## Quick Start

```julia
using CosmoDTFE

points = [Point3(rand(), rand(), rand()) for _ in 1:50]
weights = rand(50)

estimator = DensityEstimator(points, weights; depth=9)
densityValue = estimator(Point3(0.5, 0.5, 0.5))

xs = range(0.0, 1.0, length=16)
ys = range(0.0, 1.0, length=16)
zs = range(0.0, 1.0, length=16)
densityGrid = estimator((xs, ys, zs))
```

## Core API

- `Point3`: static 3D point alias, `SVector{3,Float64}`
- `DensityEstimator`: scalar DTFE density estimation
- `VelocityEstimator`: vector velocity estimation and `velocityGradient`
- `tessellate`: TetGen-backed tetrahedralization
- `BoundingVolumeHierarchy`: BVH acceleration structure

## Development Notes

The package intentionally does not track `Manifest.toml`; applications and
analysis projects should manage their own manifests. Heavy data-loading,
plotting, and Illustris workflows belong in downstream scripts or examples, not
in the package dependency surface.
