"""
    JuliaDTFE

Delaunay Tessellation Field Estimator (DTFE) implemented in Julia.

Provides high-performance density field estimation from particle data using
Delaunay tessellation with BVH-accelerated spatial queries.

# Core Features
- 3D Delaunay tessellation via TetGen
- BVH tree for fast spatial queries
- Multi-threaded field estimation

# Main Types
- `Point3`: 3D point type (SVector{3, Float64})
- `Triangulation3D`: Tessellation with density estimates
- `BoundingVolumeHierarchy`: BVH for fast lookups

# Main Functions
- `DensityEstimator`: Build DTFE estimator from point cloud
- `dtfe`: Interpolate density at a point (called by DensityEstimator functor)
"""
module JuliaDTFE

using StaticArrays
using TetGen
using Quickhull
using LinearAlgebra

# Core components (plain includes, no submodules)
include("Elements.jl")
include("BVH.jl")
include("Tessellate.jl")
include("Searchers.jl")
include("Estimators.jl")

# Public API exports
export Point3, Triangulation3D
export BoundingVolumeHierarchy, BVHTree, BVHLeaf, BVHNode
export tessellate, tessellateQH
export findSimplex, findId, recursiveSearch
export dtfe, dtfeMultiThread, DensityEstimator
export VelocityEstimator, velocityGradient

end # module
