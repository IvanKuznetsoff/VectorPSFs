module VectorPsfs

using HCubature
using StaticArrays
using LinearAlgebra
using Statistics
using SmoothingSplines
using Optim

# Include submodules (snake_case filenames)
include("nvspectrum.jl")
include("objectives.jl")
include("plane_parallel_plates.jl")
include("aberration_functions.jl")
include("psf_core.jl")
include("nvcenter.jl")
include("strehl.jl")

# Re-export or export only the public-facing API
export Objective, PlaneParallelPlate, NVCenter
export MplanApo100x, Lmplfln100xbd, MplanApo50x, M10x
export Diamond, FusedSilica
export Psf, Strehl

end # module VectorPsfs
