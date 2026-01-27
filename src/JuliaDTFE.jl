"""
    JuliaDTFE

Delaunay Tessellation Field Estimator (DTFE) implemented in Julia.

Provides high-performance density field estimation from particle data using
Delaunay tessellation with BVH-accelerated spatial queries.

# Core Features
- 3D Delaunay tessellation via TetGen
- BVH tree for fast spatial queries
- Multi-threaded field estimation
- GPU acceleration support (experimental)

# Main Types
- `Point3`: 3D point type (SVector{3, Float64})
- `Tetrahedron`: Tetrahedron with vertices and volume
- `Triangulation3D`: Tessellation with density estimates
- `BoundingVolumeHierarchy`: BVH for fast lookups

# Main Functions
- `standardEstimator`: Build DTFE estimator from point cloud
- `dtfe`: Interpolate density at a point
- `dtfeMultiThread`: Multi-threaded field estimation on grid
"""
module JuliaDTFE

using StaticArrays
using TetGen
using LinearAlgebra
using KernelAbstractions
using CUDA

# Core components (plain includes, no submodules)
include("Elements.jl")
include("BVH.jl")
include("Tessellate.jl")
include("Searchers.jl")
include("Estimators.jl")

# Public API exports
export Point3, Tetrahedron, Triangulation3D
export BoundingVolumeHierarchy, BVHTree, BVHLeaf, BVHNode
export tessellate
export findSimplex, findId, recursiveSearch
export dtfe, dtfeMultiThread, DensityEstimator

end # module
