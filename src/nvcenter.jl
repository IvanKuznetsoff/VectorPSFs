# nvcenter.jl
# Provides an `NVCenter` type to store an NV emission spectrum with associated spectral weights,
# plus methods to compute a polychromatic PSF from NV emission.

"""
    NVCenter

A struct for storing an NV emission spectrum, specifying:

- `λs::AbstractVector{<:Real}`: the sampling wavelengths in micrometers (µm)
- `weight::AbstractVector{<:Real}`: weighting (normalized spectral intensity)

Example usage:
```julia
# Construct from an existing array of wavelengths & weights
center = NVCenter([0.632, 0.64, 0.65], [0.3, 0.5, 0.2])

# Or use the fitted NV default (pulling from `nv_spectrum_spline`):
center2 = NVCenter(0.5:0.01:0.7)
```
"""
struct NVCenter
    λs::AbstractVector{<:Real}      # sampling wavelengths (µm)
    weight::AbstractVector{<:Real}  # weighting on spectral intensity
end

"""
    NVCenter(λs::Real, weight::Real)

Constructor that takes arrays of wavelengths `λs` and weights `weight` (must be same length).
Converts them to `Float64` and returns a new `NVCenter`.
"""
function NVCenter(λs::Real, weight::Real)
    @assert length(λs) == length(weight) "λs and weight must have the same length."
    return NVCenter(Float64.(λs), Float64.(weight))
end

"""
    NVCenter(λs::AbstractVector{<:Real})

Constructor that uses the global NV spectrum spline (`nv_spectrum_spline`) 
to obtain weights for the given array of wavelengths `λs`.
Weights are normalized (so their sum is 1).
"""
function NVCenter(λs::AbstractVector{<:Real})
    w = SmoothingSplines.predict(NVspectrum.nv_spectrum_spline, λs)
    weight = w ./ sum(w)
    return NVCenter(λs, weight)
end

# -------------------------------------------------------------------
# Weighted NV-luminescence PSF
# -------------------------------------------------------------------

"""
    PSF(x, y, z, obj::Objective, nv::NVCenter; rtol=1e-4, atol=1e-5)

Compute a polychromatic PSF for an NV center (no plate scenario).  
Integrates over all wavelengths in `nv.λs`, calling the single-wavelength `PSF(...)` for each λ
and weighting by `nv.weight`.

- `(x, y, z)` in µm
- `obj` is an `Objective`
- `nv::NVCenter` is the NV spectrum data
- `rtol, atol` are integration tolerances
"""
function PSF(x, y, z, obj::Objective, nv::NVCenter; rtol=1e-4, atol=1e-5)
    # Example: λ in µm
    res = map(λ -> PSF(x, y, z, λ, obj; rtol=rtol, atol=atol), nv.λs)
    return dot(res, nv.weight)
end

"""
    PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate; rtol=1e-4, atol=1e-5)

Compute a polychromatic PSF for NV emission with a normal-incidence plate.  
Broadcasts over `nv.λs`, calling single-wavelength `PSF(...)` for each λ, then weights with `nv.weight`.
"""
function PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate; rtol=1e-4, atol=1e-5)
    res = map(λ -> PSF(x, y, z, λ, obj, plate; rtol=rtol, atol=atol), nv.λs)
    return dot(res, nv.weight)
end

"""
    PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate, α::Float64; rtol=1e-4, atol=1e-5)

Compute a polychromatic PSF for NV emission with a tilted-incidence plate.  
Again, integrates over `nv.λs` and sums with `nv.weight`.
"""
function PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate, α::Float64; rtol=1e-4, atol=1e-5)
    res = map(λ -> PSF(x, y, z, λ, obj, plate, α; rtol=rtol, atol=atol), nv.λs)
    return dot(res, nv.weight)
end
