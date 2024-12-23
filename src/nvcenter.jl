# nvcenter.jl

# Struct NVCenter stores an NV emission spectrum (possibly as a smoothing spline).
# By default, we build it from globally defined wavelength_data & spectrum_data.
import ..NVspectrum

struct NVCenter
    λs::Vector{Float64}            # sampling wavelengths
    spline::SmoothingSplines.SmoothingSpline
end


# Constructor that fits a smoothing spline to the given data.

function NVCenter(λs::Vector{Float64}, spec::Vector{Float64}, smooth_param::Float64 = 250.0)
    spl = SmoothingSplines.fit(SmoothingSplines.SmoothingSpline, λs, spec, smooth_param)
    return NVCenter(λs, spl)
end

# A default NVCenter using your global data (if desired):
const default_nvcenter = NVCenter(NVspectrum.wavelength, NVspectrum.spectrum)  # example

# -------------------------------------------------------------------
# Weighted NV-luminescence PSF
# -------------------------------------------------------------------

# PSF(x, y, z, obj::Objective, nv::NVCenter; rtol, atol)

# No plate scenario, integrates over multiple wavelengths from the NV spectrum.

function PSF(x, y, z, obj::Objective, nv::NVCenter; rtol=1e-3, atol=1e-4)
    # 1) Predict intensity weights from the spline
    w = SmoothingSplines.predict(nv.spline, nv.λs)
    weight = w ./ sum(w)
    # 2) Evaluate the “aberration-free” PSF for each λ in nv.λs/1000 if needed
    #    (or store them as nm vs. µm, your call).
    #    Example assumes λ was in nm in spline, so we convert back to (nm→m) or keep consistent:
    λarray = nv.λs .* 1e-3   # if your 'λs' is in nm, this is in μm now
    res = map(λ -> PSF(x, y, z, λ, obj, nothing; rtol=rtol, atol=atol), λarray)
    return dot(res, weight)
end


# Normal-incidence plate scenario.

function PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate; rtol=1e-3, atol=1e-4)
    w = SmoothingSplines.predict(nv.spline, nv.λs)
    weight = w ./ sum(w)
    λarray = nv.λs .* 1e-3
    res = map(λ -> PSF(x, y, z, λ, obj, plate; rtol=rtol, atol=atol), λarray)
    return dot(res, weight)
end

# Tilted-incidence scenario. Example uses your paraxial shift logic.

function PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate, α::Float64; rtol=1e-3, atol=1e-4)
    # Weighted sums
    w = SmoothingSplines.predict(nv.spline, nv.λs)
    weight = w ./ sum(w)
    λarray = nv.λs .* 1e-3

    res = map(λ -> PSF(x, y, z, λ, obj, α, plate; rtol=rtol, atol=atol), λarray)
    return dot(res, weight)
end
