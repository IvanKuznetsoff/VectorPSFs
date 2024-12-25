# aberration_functions.jl
# Aberration phase functions for normal or tilted incidence.

"""
    Φ(sρ, n)

Computes the aberration phase for **normal incidence**.

- `sρ`: the product of normalized aperture `s` and radial coordinate `ρ` in the pupil plane.
- `n`: refractive index (e.g., 1.0 for air or the immersion index).
Returns a dimensionless phase shift due to the optical path difference.
"""
function Φ(sρ, n)
    @fastmath begin
        R_2 = sρ^2
        Φv = 1/n * (1 - sqrt(1 - R_2)) - n * (1 - sqrt(1 - R_2 / n^2))
    end
    return Φv
end

"""
    Φ(sρ, ϕ, n, α)

Computes the aberration phase for **tilted incidence**.

- `sρ`: again, `s * ρ`.
- `ϕ`: azimuthal angle in the pupil plane.
- `n`: refractive index.
- `α`: tilt angle (radians).

Returns a phase shift accounting for tilt plus the normal-incidence phase.
"""
function Φ(sρ, ϕ, n, α)
    @fastmath begin
        cosα = cos(α)
        sinα = sin(α)
        cosϕ = cos(ϕ)
        Φv = -cosα / sqrt(n^2 - sinα^2) * (
                sinα * cosϕ * sρ + cosα * sqrt(1 - sρ^2)
             ) +
              (1 - n^2) / sqrt(n^2 - sinα^2) +
              sqrt(
                  n^2 - sinα^2 -
                  cosα^2 * sρ^2 +
                  sinα^2 * sρ^2 * cosϕ^2 +
                  2sinα * cosα * sρ * cosϕ * sqrt(1 - sρ^2)
              )
    end
    return Φv
end