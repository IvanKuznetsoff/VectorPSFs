# strehl.jl
# Strehl ratio computations and thickness optimization for S ≈ 0.8.

"""
    strehl(t, λ, obj::Objective, plate::Function=Diamond; 
           zrange=[-2.0, 4.0], rtol=1e-10, atol=1e-10)

Compute the Strehl ratio when a plate of thickness `t` (µm) is placed in the optical path.
- `λ` is wavelength in µm
- `obj` is the objective (e.g. `MPlanApo50x()`)
- `plate` is a constructor function for a `PlaneParallelPlate` (defaults to `Diamond`).
- `zrange` in µm is the defocus search range: the function does a 1D optimization to find the best focus.
- `rtol`, `atol` are integration tolerances for the PSF.

Returns the maximum on-axis intensity (aberrated) normalized by the unaberrated case.
Example:
```julia
S = strehl(50.0, 0.7, MPlanApo100x())  # diamond thickness=50 µm, λ=0.7 µm
```
"""
function strehl(
    t::Real,
    λ::Real,
    obj::Objective,
    plate::Function=Diamond;
    zrange::AbstractVector{<:Real}=[-2.0, 4.0],
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)
    # Unaberrated on-axis intensity at best focus
    I_0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)

    # Aberrated intensity as function of defocus z
    objective_func(z) = PSF(0.0, 0.0, z, λ, obj, plate(t); rtol=rtol, atol=atol) / I_0

    # Maximize objective_func => minimize negative
    res = optimize(z -> -objective_func(z), zrange[1], zrange[2], Brent())
    return -minimum(res)
end

# -------------------------------------------------------------------
# Tilted Incidence
# -------------------------------------------------------------------

"""
    strehl(t, λ, obj::Objective, plate::Function, α::Real; 
           zrange=[-2.0, 4.0], rtol=1e-10, atol=1e-10)

Compute the Strehl ratio for a plate of thickness `t` (µm), wavelength `λ` (µm), 
Objective `obj`, and tilt angle `α` (radians).  
The search for best defocus is done in `zrange` (µm).  
"""
function strehl(
    t::Real,
    λ::Real,
    obj::Objective,
    plate::Function,
    α::Real;
    zrange::AbstractVector{<:Real}=[-2.0, 4.0],
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)
    I_0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)
    objective_func(z) = PSF(0.0, 0.0, z, λ, obj, plate(t), α; rtol=rtol, atol=atol) / I_0

    res = optimize(z -> -objective_func(z), zrange[1], zrange[2])
    return -minimum(res)
end

# -------------------------------------------------------------------
# Fixed Defocus
# -------------------------------------------------------------------

"""
    strehl(t, λ, obj::Objective, z_est::Real, plate::Function=Diamond; rtol=1e-10, atol=1e-10)

Compute the on-axis Strehl ratio for a plate thickness `t` (µm) at a specified defocus `z_est` (µm).
No tilt assumed. 
"""
function strehl(
    t::Real,
    λ::Real,
    obj::Objective,
    z_est::Real,
    plate::Function=Diamond;
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)::Float64
    I0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)
    I_tz = PSF(0.0, 0.0, z_est, λ, obj, plate(t); rtol=rtol, atol=atol)
    return I_tz / I0
end

"""
    strehl(t, λ, obj::Objective, z_est::Real, plate::Function, α::Real; rtol=1e-10, atol=1e-10)

Same as above, but with incidence angle `α` (tilt).
"""
function strehl(
    t::Real,
    λ::Real,
    obj::Objective,
    z_est::Real,
    plate::Function,
    α::Real;
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)::Float64
    I0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)
    I_tz = PSF(0.0, 0.0, z_est, λ, obj, plate(t), α; rtol=rtol, atol=atol)
    return I_tz / I0
end

# -------------------------------------------------------------------
# maxtol_thick
# -------------------------------------------------------------------

"""
    maxtol_thick(λ::Real, obj::Objective; plate=Diamond, t_range=(0.0, 500.0), zrange=[-2.0,4.0], ...)

Find thickness `t` in `t_range` that yields a Strehl ratio near 0.8.  
No tilt is assumed.  
We define:
```
cost(t) = (strehl(t, λ, obj, plate; zrange=zrange) - 0.8)^2
```
and minimize with a 1D search.  
Returns the thickness that best satisfies `S ≈ 0.8`.
"""
function maxtol_thick(
    λ::Real, 
    obj::Objective,
    plate::Function=Diamond; 
    t_range::Tuple{Real,Real}=(0.0, 500.), 
    zrange::AbstractVector{<:Real}=[-2.0, 4.0], 
    rtol::Float64=1e-10, 
    atol::Float64=1e-10
)
    function cost(t)
        s_val = strehl(t, λ, obj, plate; zrange=zrange, rtol=rtol, atol=atol)
        return (s_val - 0.8)^2
    end

    t_min, t_max = t_range
    res = optimize(cost, t_min, t_max, Brent())
    return Optim.minimizer(res)
end

"""
    maxtol_thick(λ::Real, obj::Objective, plate::Function, α::Real; t_range=(0.0, 500.0), zrange=[-2.0,4.0], ...)

Same as above, but for tilted incidence angle `α`.  
```
cost(t) = (strehl(t, λ, obj, plate, α; zrange=zrange) - 0.8)^2
```
"""
function maxtol_thick(
    λ::Real, 
    obj::Objective,
    plate::Function,
    α::Real; 
    t_range::Tuple{Real,Real}=(0.0, 500.), 
    zrange::AbstractVector{<:Real}=[-2.0, 4.0], 
    rtol::Float64=1e-10, 
    atol::Float64=1e-10
)

    function cost(t)
        s_val = strehl(t, λ, obj, plate, α; zrange=zrange, rtol=rtol, atol=atol)
        return (s_val - 0.8)^2
    end

    t_min, t_max = t_range
    res = optimize(cost, t_min, t_max, Brent())
    return Optim.minimizer(res)
end
