# JuliaDTFE

[![Julia](https://img.shields.io/badge/julia-v1.9+-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/status-active-success.svg)]()

High-performance **Delaunay Tessellation Field Estimator (DTFE)** implemented in Julia.

## Overview

DTFE is a method for estimating continuous density fields from discrete point samples. Unlike traditional methods (e.g., kernel density estimation), DTFE uses the natural geometry of the point distribution through Delaunay tessellation. That way multi-scale features can be estimated more accurately without reliance on kernel size.

This makes it a more robust method for estimating density fields from point samples, especially for cosmological applications, where multi-scale morphology is vital. This document will outline the particular choices made in this implementation. 

For more information on the original formalism, see [van de Weygaert & Schaap, 2009](https://arxiv.org/abs/0708.1441).

This project took large inspiration from the [phase space implementation by Job Feldbrugge](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl?tab=readme-ov-file).

## Algorithm

1. **Delaunay Tessellation** - The input point cloud is triangulated in 3D using TetGen, producing a mesh of tetrahedra that uniquely partitions the convex hull. This is done with the `tessellate` function. The outputs are an array of all coordinates and an array of all tetrahedra. For memory efficiency, each tetrahedron is stored as an integer array of size N×4, where N is the number of tetrahedra and the 4 columns are the indices of the vertices. The actual coordinates of each vertex are given by the coordinates array.

Next, a triangulation is created, which precomputes the volumes of each tetrahedron and the local density at each vertex (the "star" values in the original notation).

Finally, a **Bounding Volume Hierarchy (BVH)** is built for fast point location queries. The BVH recursively splits the set of tetrahedra into axis-aligned bounding boxes. At each level, the tetrahedra are split along the longest axis of the current bounds (cycling x, y, z), ensuring a balanced tree structure. This reduces the search complexity for finding the tetrahedron containing a query point from O(n) to O(log n). 

All of this is initialized when calling the `DensityEstimator` constructor.

2. **Density Estimation** 

- For each query point, the method:
   - Traverses the BVH to find the leaf node containing the point. In an ideal scenario the destination leaf of the BVH tree will contain a single tetrahedron. Realistically, the destination leaf will contain a neighborhood of tetrahedra.

   - Linearly searches the tetrahedra within that leaf (checking barycentric coordinates) to find the exact enclosing tetrahedron using `earlyStopSearch`. Delaunty Tesselations guarantee that there will be only one valid tetrahedron containing the point, so we on average half operations when stopping early.

   - Computes the local density `ρ(x)` via barycentric interpolation of the densities at the four vertices of the enclosing tetrahedron:
   ```
   ρ(x) = Σ λᵢ ρ(vᵢ)
   ```
   where `λᵢ` are the barycentric coordinates and `ρ(vᵢ)` are the densities at the vertices. This is more stable and faster than gradient-based methods.

To interpotale, we simply call the density estimator functor, as seen in the quick start example.


### Key Features

- **BVH-accelerated spatial queries** - O(log n) point location for efficient lookups in large datasets.
- **Multi-threaded grid evaluation** - Parallel computation when estimating density fields on grids or large lists of points.
- **GPU acceleration** - Experimental CUDA support via KernelAbstractions.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ispirovjr/JuliaDTFE")
```

Or for local development:
```julia
] activate .
] instantiate
```

## Quick Start

```julia

using JuliaDTFE

points = [Point3(rand(), rand(), rand()) for _ in 1:50]
weights = rand(50)

est = DensityEstimator(points, weights,9)

density = est([0.5, 0.5, 0.5]) # for single point

xs = range(0.0, 1.0, length=3)
ys = range(0.0, 1.0, length=3)
zs = range(0.0, 1.0, length=3)

densityGrid = est((xs, ys, zs)) # for grid

```

## Project Structure

```
JuliaDTFE/
├── src/
│   ├── JuliaDTFE.jl      # Main module
│   ├── Elements.jl       # Core types: Point3, Tetrahedron, Triangulation3D
│   ├── BVH.jl            # Bounding Volume Hierarchy for spatial queries
│   ├── Tessellate.jl     # TetGen wrapper for Delaunay tessellation
│   ├── Searchers.jl      # Point-in-tetrahedron search algorithms
│   ├── Estimators.jl     # DTFE density estimation functions
│   └── Plotting.jl       # Visualization utilities (GLMakie)
├── test/                 # Unit tests for each module
├── examples/             # Benchmarks and usage examples
└── saves/                # Saved data files
```
