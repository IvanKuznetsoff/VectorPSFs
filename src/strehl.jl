# strehl.jl

# -----------------------------------------------------------------------------
# Computes the Strehl ratio with a diamond plate of thickness `t` at wavelength `λ`, 
# using an Objective `obj` (normal incidence). 
# The Strehl ratio is the maximum on-axis intensity normalized by the unaberrated case.
# -----------------------------------------------------------------------------
function Strehl(
    t::Float64,
    λ::Float64,
    obj::Objective,
    plate::Function=Diamond;
    zrange::Vector{Float64}=[-2.0, 4.0],
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)
    # Unaberrated on-axis intensity at best focus
    I_0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)

    # Aberrated intensity as function of defocus z
    objective_func(z) = PSF(0.0, 0.0, z, λ, obj, plate(t); rtol=rtol, atol=atol) / I_0

    # Maximize objective_func => minimize -objective_func
    res = optimize(z -> -objective_func(z), zrange[1], zrange[2], Brent())
    return -minimum(res)
end

# -----------------------------------------------------------------------------
# Computes the Strehl ratio for a diamond plate of thickness `t`, wavelength `λ`, 
# and an Objective `obj` at incidence angle `α`.
# The Strehl ratio at the best defocus in [zrange[1], zrange[2]].
# -----------------------------------------------------------------------------
function Strehl(
    t::Float64,
    λ::Float64,
    α::Float64,
    obj::Objective,
    plate::Function=Diamond;
    zrange::Vector{Float64}=[-2.0, 4.0],
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)
    # Unaberrated on-axis intensity 
    I_0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)

    # Aberrated intensity vs. defocus z
    objective_func(z) = PSF(0.0, 0.0, z, λ, obj, α, plate(t); rtol=rtol, atol=atol) / I_0

    res = optimize(z -> -objective_func(z), zrange[1], zrange[2])
    return -minimum(res)
end

# -----------------------------------------------------------------------------
# Computes the on-axis Strehl ratio for a diamond plate of thickness `t` at a 
# user-specified defocus `z_est`.
# -----------------------------------------------------------------------------
function Strehl(
    t::Float64,
    λ::Float64,
    obj::Objective,
    z_est::Float64,
    plate::Function=Diamond;
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)::Float64
    # On-axis intensity with no plate
    I0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)

    # PSF with plate(t) at that defocus
    I_tz = PSF(0.0, 0.0, z_est, λ, obj, plate(t); rtol=rtol, atol=atol)

    return I_tz / I0
end

# -----------------------------------------------------------------------------
# Same as above but with incidence angle α.
# -----------------------------------------------------------------------------
function Strehl(
    t::Float64,
    λ::Float64,
    α::Float64,
    obj::Objective,
    z_est::Float64,
    plate::Function=Diamond;
    rtol::Float64=1e-10,
    atol::Float64=1e-10
)::Float64
    # On-axis intensity with no plate
    I0 = PSF(0.0, 0.0, 0.0, λ, obj, plate(1e-3); rtol=rtol, atol=atol)

    # PSF with plate(t) at that defocus
    I_tz = PSF(0.0, 0.0, z_est, λ, obj, α, plate(t); rtol=rtol, atol=atol)

    return I_tz / I0
end

# -----------------------------------------------------------------------------
# Finds the thickness t_star such that the Strehl ratio is as close to 0.8 
# as possible, given a user-specified z_spline(t).
# -----------------------------------------------------------------------------
function maxtol_thick(
    λ::Float64, 
    obj::Objective,
    plate::Function=Diamond; 
    t_range::Tuple{Real, Real}=(0.0, 500.), 
    zrange::Vector{Float64}=[-2.0, 4.0], 
    rtol::Float64=1e-10, 
    atol::Float64=1e-10
)
    # cost(t) = [S(t) - 0.8]^2 using the defocus from z_spline
    function cost(t)
        s_val = Strehl(t, λ, obj, plate; zrange=zrange, rtol=rtol, atol=atol)
        return (s_val - 0.8)^2
    end

    t_min, t_max = t_range
    res = optimize(cost, t_min, t_max, Brent())
    return Optim.minimizer(res)
end

# -----------------------------------------------------------------------------
# Same as above but includes incidence angle α.
# -----------------------------------------------------------------------------
function maxtol_thick(
    λ::Float64, 
    obj::Objective,
    α::Float64,
    plate::Function=Diamond; 
    t_range::Tuple{Real, Real}=(0.0, 500.), 
    zrange::Vector{Float64}=[-2.0, 4.0], 
    rtol::Float64=1e-10, 
    atol::Float64=1e-10
)

    # cost(t) = [S(t) - 0.8]^2 using the defocus from z_spline
    function cost(t)
        s_val = Strehl(t, λ, α, obj, plate; zrange=zrange, rtol=rtol, atol=atol)
        return (s_val - 0.8)^2
    end

    t_min, t_max = t_range
    res = optimize(cost, t_min, t_max, Brent())
    return Optim.minimizer(res)
end