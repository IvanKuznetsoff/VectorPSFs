module VectorPSFs

using HCubature
using StaticArrays
using LinearAlgebra
using Statistics
using SmoothingSplines
using SmoothingSplines: SmoothingSpline, predict

using Optim

# Include submodules (snake_case filenames)
include("NVspectrum.jl")
include("objectives.jl")
include("plane_parallel_plates.jl")
include("aberration_functions.jl")
include("psf_core.jl")
include("nvcenter.jl")
include("strehl.jl")

# From `objectives.jl`
export Objective
export MPlanApo100x, LMPLFLN100XBD, MPlanApo50x, M10x, UPLXAPO100X

# From `plane_parallel_plates.jl`
export PlaneParallelPlate
export Diamond, FusedSilica, BorosilicateCrown, Sapphire, MagnesiumFluoride
export CustomPlate  # for user-defined Sellmeier plates

# From `nvcenter.jl`
export NVCenter  # optionally export default_nvcenter if you want

# From `psf_core.jl`
export PSF

# From `strehl.jl`
export Strehl, maxtol_thick


end # module VectorPsfs
