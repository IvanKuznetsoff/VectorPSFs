# plane_parallel_plates.jl

"""
PlaneParallelPlate stores a plate thickness `t` (in nanometers)
and a function `n_λ` giving the refractive index at wavelength λ (in micrometers).
"""
struct PlaneParallelPlate
    t::Float64        # Thickness in nanometers
    n_λ::Function     # Refractive Index Curve, n(λ in µm)
end


# ------------------------------------------------------------
# Diamond
# ------------------------------------------------------------
"""
Diamond(t)

Approximate dispersion for diamond.
`λ` is assumed in micrometers (µm).
Internally converts µm → nm by multiplying by 1e3.
"""
function Diamond(t)
    n_λ(λ) = sqrt(1 + (4.658 * (λ * 1e3)^2) /
                         ((λ * 1e3)^2 - 112.5^2))
    return PlaneParallelPlate(t, n_λ)
end


# ------------------------------------------------------------
# Sellmeier Utility
# ------------------------------------------------------------
"""
sellmeier(λ, B, C)

Generic Sellmeier equation:
    n^2 = 1 + ∑[ Bᵢ * λ² / (λ² - Cᵢ) ]
where `λ` is the wavelength in micrometers (µm),
and `B`/`C` are vectors of Sellmeier coefficients.
"""
function sellmeier(λ, B::Vector{Float64}, C::Vector{Float64})
    @assert length(B) == length(C) "B and C must have the same length."
    λ2 = λ^2
    val = 0.0
    @inbounds @simd for i in eachindex(B)
        val += B[i] * λ2 / (λ2 - C[i])
    end
    return sqrt(1 + val)
end


# ------------------------------------------------------------
# Fused Silica (SiO₂): Malitson
#   Source: https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson
# ------------------------------------------------------------
function FusedSilica(t)
    B = [0.6961663, 0.4079426, 0.8974794]
    C = [0.0684043^2, 0.1162414^2, 9.896161^2]

    n_λ(λ) = sellmeier(λ, B, C)
    return PlaneParallelPlate(t, n_λ)
end


# ------------------------------------------------------------
# Borosilicate Crown (BK7): Schott
#   Source: https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
# ------------------------------------------------------------
function BorosilicateCrown(t)
    # Schott BK7 Sellmeier coefficients
    B = [1.040, 0.230, 1.010]
    C = [0.006, 0.020, 103.560]

    n_λ(λ) = sellmeier(λ, B, C)
    return PlaneParallelPlate(t, n_λ)
end


# ------------------------------------------------------------
# Sapphire (Al₂O₃) ordinary ray: Li
#   Source: https://refractiveindex.info/?shelf=main&book=Al2O3&page=Li
# ------------------------------------------------------------
function Sapphire(t)
    # Ordinary-ray Sellmeier coefficients
    B = [1.4313493, 0.65054713, 5.3414021]
    C = [0.0726631^2, 0.1193242^2, 18.028251^2]

    n_λ(λ) = sellmeier(λ, B, C)
    return PlaneParallelPlate(t, n_λ)
end


# ------------------------------------------------------------
# Magnesium Fluoride (MgF₂) ordinary ray: Li
#   Source: https://refractiveindex.info/?shelf=main&book=MgF2&page=Li
# ------------------------------------------------------------
function MagnesiumFluoride(t)
    # Ordinary-ray Sellmeier coefficients
    B = [0.48755108, 0.39875031, 2.3120353]
    C = [0.04338408^2, 0.09461442^2, 23.793604^2]

    n_λ(λ) = sellmeier(λ, B, C)
    return PlaneParallelPlate(t, n_λ)
end

# ------------------------------------------------------------
# Custom Plate
# ------------------------------------------------------------
function CustomPlate(t::Float64, B::Vector{Float64}, C::Vector{Float64}) 
    n_λ(λ) = sellmeier(λ, B, C) 
    return PlaneParallelPlate(t, n_λ) 
end