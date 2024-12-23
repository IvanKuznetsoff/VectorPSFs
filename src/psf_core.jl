# psf_core.jl

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
function PSF(x, y, z, λ, NA; rtol=1e-3, atol=1e-4)
    n = 1.0
    s = NA
    k = 2π / λ
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)
    return psf_inner(ψ, u, v, s, n; rtol=rtol, atol=atol)
end

# -------------------------------------------------------------------
# Aberration Free + Objective (no plate => aberration free w/ immersion)
# -------------------------------------------------------------------
function PSF(x, y, z, λ, obj::Objective;
    rtol::Float64=1e-3, atol::Float64=1e-4)
    # If there's no plate, we might rely on simpler code that doesn't 
    # set n_a, t, etc. E.g.:

    n_imm = obj.n
    s     = obj.NA / obj.n
    k     = 2π * n_imm / λ

    # For "normal incidence" on the sample side (no plate),
    # let n = 1 (air or sample index if needed).
    # The defocus factor:
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2)) * z
    #  ^ note that  s^2 / n^2 = s^2 / 1^2 = s^2

    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)

    # Then call the "aberration-free" psf_inner(...) 
    # from your psf_core logic:
    return psf_inner(ψ, u, v, s, 1.0; rtol=rtol, atol=atol)
end

# -------------------------------------------------------------------
# Normal incidence. Keep the old logic (air, n=1).
# -------------------------------------------------------------------

function PSF(x, y, z, λ, NA::Float64, plate::PlaneParallelPlate; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
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

# -------------------------------------------------------------------
# Normal incidence + Objective
# -------------------------------------------------------------------

function PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
    n_imm = obj.n
    s     = obj.NA / obj.n
    k     = 2π * n_imm / λ       # n factor in immersion
    n     = target.n_λ(λ)        # plate’s n(λ), typically
    n_a   = plate.n_λ(λ) / n_imm # n divided by n_imm for Φ
    t     = plate.t
    # Same integration logic
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)
    return psf_inner(ψ, u, v, k, s, n, n_a, t; rtol=rtol, atol=atol)
end

# -------------------------------------------------------------------
# Tilted incidence. Keep the old logic (air, n=1).
# -------------------------------------------------------------------

function PSF(x, y, z, λ, NA::Float64, plate::PlaneParallelPlate, α; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
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

# -------------------------------------------------------------------
# Tilted incidence + Objective
# -------------------------------------------------------------------

function PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate, α; target::PlaneParallelPlate=plate, rtol=1e-3, atol=1e-4)
    n_imm = obj.n
    s     = obj.NA / obj.n
    k     = 2π * n_imm / λ
    n     = target.n_λ(λ)
    n_a   = plate.n_λ(λ) / n_imm  # (2) 
    t     = plate.t
    u = k * (1 - sqrt(1 - s^2)) / (1 - sqrt(1 - s^2 / n^2)) * z
    v = k * s * sqrt(x^2 + y^2)
    ψ = atan(x, y)
    return psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol=rtol, atol=atol)
end


abstract type PSFMode end

struct NoAberration <: PSFMode end
struct NormalIncidence <: PSFMode end
struct TiltedIncidence <: PSFMode end


mutable struct PSFParams{M<:PSFMode}
    λ::Float64                    # wavelength in µm
    obj::Objective               # objective lens
    plate::Union{PlaneParallelPlate,Nothing}  # nothing if no plate
    α::Float64                   # tilt angle (radians), 0.0 if no tilt
end

# -------------------------------------------------------------------
# Aberration Free + Objective (no plate => aberration free w/ immersion)
# -------------------------------------------------------------------

function PSFParams(λ::Float64, obj::Objective)
    return PSFParams{NoAberration}(λ, obj, nothing, 0.0)
end

function PSF(x, y, z, 
    param::PSFParams{NoAberration};
    rtol::Float64=1e-3, atol::Float64=1e-4)

    return PSF(x, y, z, param.λ, param.obj; rtol=rtol, atol=atol)
end

# -------------------------------------------------------------------
# Normal incidence, with plane-parallel plate
# -------------------------------------------------------------------
function PSFParams(λ::Float64, obj::Objective, plate::PlaneParallelPlate)
    return PSFParams{NormalIncidence}(λ, obj, plate, 0.0)
end

function PSF(x, y, z, 
    param::PSFParams{NormalIncidence};
    rtol::Float64=1e-3, atol::Float64=1e-4)

    return PSF(x, y, z, param.λ, param.obj, param.plate; rtol=rtol, atol=atol)
end

# -------------------------------------------------------------------
# Tilted incidence, with plane-parallel plate
# -------------------------------------------------------------------
function PSFParams(λ::Float64, obj::Objective, plate::PlaneParallelPlate, α::Float64)
    return PSFParams{TiltedIncidence}(λ, obj, plate, α)
end

function PSF(x, y, z, 
    param::PSFParams{TiltedIncidence};
    rtol::Float64=1e-3, atol::Float64=1e-4)

    return PSF(x, y, z, param.λ, param.obj, param.plate, param.α; rtol=rtol, atol=atol)
end
