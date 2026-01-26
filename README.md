# JuliaDTFE

High-performance **Delaunay Tessellation Field Estimator (DTFE)** implemented in Julia.

## Overview

DTFE is a method for estimating continuous density fields from discrete point samples. Unlike traditional methods (e.g., kernel density estimation), DTFE uses the natural geometry of the point distribution through Delaunay tessellation.

### Algorithm

1. **Delaunay Tessellation** - The input point cloud is triangulated in 3D using TetGen, producing a mesh of tetrahedra that uniquely partitions the convex hull.

2. **Density Estimation** - For each point, the local density ρ* is computed as:
   ```
   ρ*ᵢ = Σⱼ (wᵢ / Vⱼ)
   ```
   where the sum is over all tetrahedra containing point i, wᵢ is the weight (mass) at that point, and Vⱼ is the tetrahedron volume.

3. **Interpolation** - To query density at any location, we find the enclosing tetrahedron using a BVH tree, then use barycentric interpolation of the vertex densities.

### Key Features

- **BVH-accelerated spatial queries** - O(log n) point location instead of O(n)
- **Multi-threaded grid evaluation** - Parallel density field computation
- **GPU acceleration** - Experimental CUDA support via KernelAbstractions

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/your-username/JuliaDTFE")
```

Or for local development:
```julia
] activate .
] instantiate
```

## Quick Start

```julia
using JuliaDTFE
using StaticArrays

# Create random point cloud
points = [Point3(rand(), rand(), rand()) for _ in 1:1000]

# Build DTFE estimator
bvh, triangulation, tetrahedra = standardEstimator(points)

# Query density at a point
density = dtfe([0.5, 0.5, 0.5], bvh, tetrahedra, triangulation)

# Evaluate on a 3D grid (multi-threaded)
xs = range(0.0, 1.0, length=50)
ys = range(0.0, 1.0, length=50)
zs = range(0.0, 1.0, length=50)
densityField = dtfeMultiThread((xs, ys, zs), bvh, tetrahedra, triangulation)
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

## API Reference

### Types

| Type | Description |
|------|-------------|
| `Point3` | 3D point (`SVector{3, Float64}`) |
| `Tetrahedron` | Tetrahedron with vertices and precomputed volume |
| `Triangulation3D` | Tessellation result with per-vertex densities |
| `BoundingVolumeHierarchy` | Spatial index for fast lookups |

### Functions

| Function | Description |
|----------|-------------|
| `tessellate(points)` | Compute Delaunay tessellation |
| `standardEstimator(points, [weights], [depth])` | Build DTFE estimator |
| `dtfe(point, bvh, tets, triangulation)` | Query density at single point |
| `dtfeMultiThread(grid, bvh, tets, triangulation)` | Parallel grid evaluation |
| `findId(point, simplices, bvh)` | Find tetrahedron containing point |

## Performance

- BVH depth controls query speed vs. memory tradeoff (default: 9)
- Multi-threading scales with available cores (`Threads.nthreads()`)
- Uses bitwise `&` in intersection tests to avoid branch prediction overhead

## License

MIT
