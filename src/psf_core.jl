# psf_core.jl

import ..phase_functions: Φ

# -------------------------------------------------------------------
# Aberration-free phase function
# -------------------------------------------------------------------
function Ξ(ρ, ϕ, ψ, u, v, s, n)
    @fastmath Ξv = v * ρ * cos(ϕ - ψ) + u * sqrt(1 - s^2 * ρ^2 / n^2)
end

function ξx(ρ, ϕ, Ξ, s)
    @fastmath begin
        R_2 = s^2 * ρ^2
        res = exp(im * Ξ) * ρ * s^2 *
              ( sqrt(1 - R_2) +
                (1 - sqrt(1 - R_2)) * (1 - cos(2ϕ)) / 2
              ) / (1 - R_2)^(1/4)
    end
    return res
end

function ξy(ρ, ϕ, Ξ, s)
    @fastmath begin
        R_2 = s^2 * ρ^2
        res = exp(im * Ξ) * ρ * s^2 *
              (1 - sqrt(1 - R_2)) * sin(2ϕ) / 2 / (1 - R_2)^(1/4)
    end
    return res
end

function ξz(ρ, ϕ, Ξ, s)
    @fastmath begin
        R_2 = s^2 * ρ^2
        res = exp(im * Ξ) * R_2 * s * cos(ϕ) / 2 / (1 - R_2)^(1/4)
    end
    return res
end

# Inner integral for the aberration-free PSF
function psf_inner(ψ, u, v, s, n; rtol=1e-3, atol=1e-4)
    resx = hcubature(
        dr -> ξx(dr[1], dr[2], Ξ(dr[1], dr[2], ψ, u, v, s, n), s),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )
    resy = hcubature(
        dr -> ξy(dr[1], dr[2], Ξ(dr[1], dr[2], ψ, u, v, s, n), s),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )
    resz = hcubature(
        dr -> ξz(dr[1], dr[2], Ξ(dr[1], dr[2], ψ, u, v, s, n), s),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )

    ret = 0.0
    ret += abs2(resx[1]::ComplexF64)
    ret += abs2(resy[1]::ComplexF64)
    ret += abs2(resz[1]::ComplexF64)
    return ret / π^2
end

# -------------------------------------------------------------------
# Normal incidence with a plane-parallel plate
# -------------------------------------------------------------------
function Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
    Φv = t * Φ(s * ρ, n_a)
    @fastmath Ξv = k * Φv + Ξ(ρ, ϕ, ψ, u, v, s, n)
end

function ξx(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
    Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
    return ξx(ρ, ϕ, Ξv, s)
end

function ξy(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
    Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
    return ξy(ρ, ϕ, Ξv, s)
end

function ξz(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
    Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
    return ξz(ρ, ϕ, Ξv, s)
end

function psf_inner(ψ, u, v, k, s, n, n_a, t; rtol=1e-3, atol=1e-4)
    resx = hcubature(
        dr -> ξx(dr[1], dr[2], ψ, u, v, k, s, n, n_a, t),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )
    resy = hcubature(
        dr -> ξy(dr[1], dr[2], ψ, u, v, k, s, n, n_a, t),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )
    resz = hcubature(
        dr -> ξz(dr[1], dr[2], ψ, u, v, k, s, n, n_a, t),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )

    ret = 0.0
    ret += abs2(resx[1]::ComplexF64)
    ret += abs2(resy[1]::ComplexF64)
    ret += abs2(resz[1]::ComplexF64)
    return ret / π^2
end

# -------------------------------------------------------------------
# Tilted incidence with a plane-parallel plate
# -------------------------------------------------------------------
function Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
    Φv = t * Φ(s * ρ, ϕ, n_a, α)
    @fastmath Ξv = k * Φv + Ξ(ρ, ϕ, ψ, u, v, s, n)
end

function ξx(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
    Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
    return ξx(ρ, ϕ, Ξv, s)
end

function ξy(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
    Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
    return ξy(ρ, ϕ, Ξv, s)
end

function ξz(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
    Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
    return ξz(ρ, ϕ, Ξv, s)
end

function psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol=1e-3, atol=1e-4)
    resx = hcubature(
        dr -> ξx(dr[1], dr[2], ψ, u, v, k, s, α, n, n_a, t),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )
    resy = hcubature(
        dr -> ξy(dr[1], dr[2], ψ, u, v, k, s, α, n, n_a, t),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )
    resz = hcubature(
        dr -> ξz(dr[1], dr[2], ψ, u, v, k, s, α, n, n_a, t),
        [0.0, 0.0], [1.0, 2π],
        rtol=rtol, atol=atol
    )

    ret = 0.0
    ret += abs2(resx[1]::ComplexF64)
    ret += abs2(resy[1]::ComplexF64)
    ret += abs2(resz[1]::ComplexF64)
    return ret / π^2
end

# -------------------------------------------------------------------
# High-level PSF(...) methods
# -------------------------------------------------------------------
# If no objective is provided, keep the old logic (air, n=1).
function PSF(x, y, z, λ, NA, plate::Nothing; target=plate, rtol=1e-3, atol=1e-4)
    n = 1.0
    s = NA
    k = 2π / λ
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)
    return psf_inner(ψ, u, v, s, n; rtol=rtol, atol=atol)
end

# Objective with no plate -> incorporate immersion n_obj
function PSF(x, y, z, λ, obj::Objective, plate::Nothing; rtol=1e-3, atol=1e-4)
    n_obj = obj.n        # immersion liquid refractive index
    s     = obj.NA
    k     = 2π * n_obj / λ  # (1) use n_obj
    # Since there's no plate, we still consider sample in air: n=1
    # The s^2/(n^2) factor becomes s^2 / 1^2
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2)) * z
    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)
    return psf_inner(ψ, u, v, s, 1.0; rtol=rtol, atol=atol)
end

# Plane-parallel plate, but no explicit Objective struct
function PSF(x, y, z, λ, NA, plate::PlaneParallelPlate; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
    n   = target.n_λ(λ)
    n_a = plate.n_λ(λ)  # remains unscaled if no objective
    t   = plate.t
    s   = NA
    k   = 2π / λ
    u   = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v   = k * s * sqrt(x^2 + y^2)
    ψ   = atan(x, y)
    return psf_inner(ψ, u, v, k, s, n, n_a, t; rtol=rtol, atol=atol)
end

# Objective with plane-parallel plate
function PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
    n_obj = obj.n
    s     = obj.NA
    k     = 2π * n_obj / λ       # n factor in immersion
    n     = target.n_λ(λ)        # plate’s n(λ), typically
    n_a   = plate.n_λ(λ) / n_obj # n divided by n_obj for Φ
    t     = plate.t
    # Same integration logic
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)
    return psf_inner(ψ, u, v, k, s, n, n_a, t; rtol=rtol, atol=atol)
end

# Tilted incidence: no explicit Objective
function PSF(x, y, z, λ, NA, α, plate::PlaneParallelPlate; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
    n   = target.n_λ(λ)
    n_a = plate.n_λ(λ)
    t   = plate.t
    s   = NA
    k   = 2π / λ
    u   = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v   = k * s * sqrt(x^2 + y^2)
    ψ   = atan(x, y)
    return psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol=rtol, atol=atol)
end

# Tilted incidence + Objective
function PSF(x, y, z, λ, obj::Objective, α, plate::PlaneParallelPlate; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
    n_obj = obj.n
    s     = obj.NA
    k     = 2π * n_obj / λ
    n     = target.n_λ(λ)
    n_a   = plate.n_λ(λ) / n_obj  # (2) 
    t     = plate.t
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)
    return psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol=rtol, atol=atol)
end