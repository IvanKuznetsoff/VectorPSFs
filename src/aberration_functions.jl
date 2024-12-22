# aberration_functions.jl

# Tilted and normal incidence phase functions (Φ)

# Normal Incidence
function Φ(sρ, n)
    @fastmath begin
        R_2 = sρ^2
        Φv = 1/n * (1 - sqrt(1 - R_2)) - n * (1 - sqrt(1 - R_2 / n^2))
    end
    return Φv
end

# Tilted Incidence
function Φ(sρ, ϕ, n, α)
    @fastmath begin
        cosα = cos(α)
        sinα = sin(α)
        cosϕ = cos(ϕ)
        Φv = -cosα / sqrt(n^2 - sinα^2) * (
                sinα * cosϕ * sρ + cosα * sqrt(1 - sρ^2)
             ) +
              (1 - n^2) / sqrt(n^2 - sinα^2) +
              sqrt(n^2 - sinα^2 -
                   cosα^2 * sρ^2 +
                   sinα^2 * sρ^2 * cosϕ^2 +
                   2sinα * cosα * sρ * cosϕ * sqrt(1 - sρ^2))
    end
    return Φv
end
