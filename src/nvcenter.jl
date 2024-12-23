# nvcenter.jl

# Struct NVCenter stores an NV emission spectrum (possibly as a smoothing spline).
# By default, we build it from globally defined wavelength_data & spectrum_data.

struct NVCenter
    λs::Vector{Float64}            # sampling wavelengths (microometer)
    weight::Vector{Float64}            # weighting on spectral intensity
end

# Constructor that fits a smoothing spline to the given data.

function NVCenter(λs::Real, weight::Real)
    @assert length(λs) == length(weight)
    return NVCenter(Float64.(λs), Float64.(weight))
end

# If only λ is specified, then give a spline-interpolated weight curve
function NVCenter(λs::Vector{Float64})
    w = SmoothingSplines.predict(nv_spectrum_spline, λs)
    weight = w ./ sum(w)
    return NVCenter(λs, weight)
end

# -------------------------------------------------------------------
# Weighted NV-luminescence PSF
# -------------------------------------------------------------------

# PSF(x, y, z, obj::Objective, nv::NVCenter; rtol, atol)

# No plate scenario, integrates over multiple wavelengths from the NV spectrum.

function PSF(x, y, z, obj::Objective, nv::NVCenter; rtol=1e-3, atol=1e-4)
    #    Example assumes λ was in nm in spline, so we convert back to (nm→m) or keep consistent:
    res = map(λ -> PSF(x, y, z, λ, obj, nothing; rtol=rtol, atol=atol), nv.λs)
    return dot(res, nv.weight)
end


# Normal-incidence plate scenario.

function PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate; rtol=1e-3, atol=1e-4)
    res = map(λ -> PSF(x, y, z, λ, obj, plate; rtol=rtol, atol=atol), nv.λs)
    return dot(res, nv.weight)
end

# Tilted-incidence scenario. 

function PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate, α::Float64; rtol=1e-3, atol=1e-4)
    res = map(λ -> PSF(x, y, z, λ, obj, α, plate; rtol=rtol, atol=atol), nv.λs)
    return dot(res, nv.weight)
end
