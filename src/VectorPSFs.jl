module VectorPSFs
# Write your package code here.
    ## load NV spectrum
    using HCubature
    using StaticArrays
    using LinearAlgebra
    using Statistics
    using SmoothingSplines
    using Optim
    include("NVSpectrum.jl")
    struct Objective
        f::Float64
        NA::Float64
        n::Float64
    end
    function MPlanApo100x()
        return Objective(2000, 0.7, 1.)
    end
    function LMPLFLN100xBD()
        return Objective(1800, 0.8, 1.)
    end
    function MPlanApo50x()
        return Objective(4000, 0.55, 1.)
    end
    function M10x()
        return Objective(16500, 0.25, 1.)
    end
    struct PlaneParallelPlate
        t::Float64 # Thickness
        n_λ::Function # Refractive Index Curve
    end
    function Diamond(t)
        n_λ(λ) = sqrt(1 + (4.658 * (λ * 1e3)^2) / ((λ * 1e3)^2 - 112.5^2))
        return PlaneParallelPlate(t, n_λ)
    end
    function FusedSilica(t)
        n_λ(λ) = sqrt(1 + (4.658 * (λ * 1e3)^2) / ((λ * 1e3)^2 - 112.5^2))
        return PlaneParallelPlate(t, n_λ)
    end
    function Φ(sρ, n)
        @fastmath begin
        R_2 =  sρ^2
        Φv = 1/n * (1 - sqrt(1 - R_2)) - n *(1 - sqrt(1-R_2/n^2))
        end
        return Φv
    end
    function Φ(sρ, ϕ, n, α)
        @fastmath begin
        cosα = cos(α)
        sinα = sin(α)
        cosϕ = cos(ϕ)
        Φv = -cosα/sqrt(n^2 - sinα^2) * (sinα * cosϕ * sρ + cosα * sqrt(1 - sρ^2)) + 
            (1 - n^2)/sqrt(n^2 - sinα^2) + 
            sqrt(
                n^2 - sinα^2 -cosα^2*sρ^2 + 
                sinα^2*sρ^2*cosϕ^2 + 
                2*sinα*cosα*sρ*cosϕ*sqrt(1-sρ^2)
            )
        end
        return Φv
    end
    function ξx(ρ, ϕ, Ξ, s)
        @fastmath begin
            R_2 =  s^2 * ρ^2
            res = exp(im*Ξ) * ρ * s^2 * (sqrt(1-R_2) + (1 - sqrt(1-R_2)) * (1 - cos(2ϕ))/2)/(1-R_2)^(1/4)
        end
        return res
    end
    function ξy(ρ, ϕ, Ξ, s)
        @fastmath begin
            R_2 =  s^2 * ρ^2
            res = exp(im*Ξ) * ρ * s^2 * (1 - sqrt(1-R_2)) * sin(2ϕ)/2/(1-R_2)^(1/4)
        end
        return res
    end
    function ξz(ρ, ϕ, Ξ, s)
        @fastmath begin
            R_2 =  s^2 * ρ^2
            res = exp(im*Ξ) * R_2 * s * cos(ϕ)/2/(1-R_2)^(1/4)
        end
        return res
    end
    # Aberration Free
    function Ξ(ρ, ϕ, ψ, u, v, s, n)
        @fastmath Ξv = v * ρ * cos(ϕ - ψ) + u * sqrt(1- s^2 * ρ^2/n^2)
    end
    function ξx(ρ, ϕ, ψ, u, v, s, n)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, s, n)
        ret = ξx(ρ, ϕ, Ξv, s)
        return ret
    end
    function ξy(ρ, ϕ, ψ, u, v, s, n)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, s, n)
        ret = ξy(ρ, ϕ, Ξv, s)
        return ret
    end
    function ξz(ρ, ϕ, ψ, u, v, s, n)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, s, n)
        ret = ξz(ρ, ϕ, Ξv, s)
        return ret
    end
    function psf_inner(ψ, u, v, s, n; rtol = 1e-3, atol = 1e-4)
        resx = hcubature(
            dr -> ξx(dr[1], dr[2], ψ, u, v, s, n),
            [0.0,0.0],
            [1.0,2π],
            rtol = rtol,atol = atol)
        resy = hcubature(
            dr -> ξy(dr[1], dr[2], ψ, u, v, s, n),
                [0.0,0.0],
                [1.0,2π],
            rtol = rtol,atol = atol)
        resz = hcubature(
            dr -> ξz(dr[1], dr[2], ψ, u, v, s, n),
                [0.0,0.0],
                [1.0,2π],
            rtol = rtol,atol = atol)
        ret = 0.0
        ret += abs2(resx[1]::ComplexF64)
        ret += abs2(resy[1]::ComplexF64)
        ret += abs2(resz[1]::ComplexF64)
        return  ret/ pi^2
    end
    # Normal Incidence
    function Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
        Φv = t * Φ(s * ρ, n_a)
        @fastmath Ξv = k * Φv + Ξ(ρ, ϕ, ψ, u, v, s, n)
    end
    function Ξ(
            ρ, ϕ, ψ, u, v, k, s, n, 
            n_a::SVector{N,T}, 
            t::SVector{N,T}) where {N, T<:Real}
        Φv = @. t * Φ(s * ρ, n_a)
        @fastmath Ξv = k * sum(Φv) + Ξ(ρ, ϕ, ψ, u, v, s, n)
    end
    function ξx(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
        ret = ξx(ρ, ϕ, Ξv, s)
        return ret
    end
    function ξy(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
        ret = ξy(ρ, ϕ, Ξv, s)
        return ret
    end
    function ξz(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)
        ret = ξz(ρ, ϕ, Ξv, s)
        return ret
    end
    function psf_inner(ψ, u, v, k, s, n, n_a, t; rtol = 1e-3, atol = 1e-4)
        resx = hcubature(
            dr -> ξx(dr[1], dr[2], ψ, u, v, k, s, n, n_a, t), 
            [0.0,0.0],
            [1.0,2π],
            rtol = rtol,atol = atol)
        resy = hcubature(
            dr -> ξy(dr[1], dr[2], ψ, u, v, k, s, n, n_a, t), 
                [0.0,0.0],
                [1.0,2π],
            rtol = rtol,atol = atol)
        resz = hcubature(
            dr -> ξz(dr[1], dr[2], ψ, u, v, k, s, n, n_a, t), 
                [0.0,0.0],
                [1.0,2π],
            rtol = rtol,atol = atol)
        ret = 0.0
        ret += abs2(resx[1]::ComplexF64)
        ret += abs2(resy[1]::ComplexF64)
        ret += abs2(resz[1]::ComplexF64)
        return  ret/ pi^2
    end
    # Tilted Incidence
    function Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
        Φv = t * Φ(s*ρ, ϕ, n_a, α)
        @fastmath Ξv = k * Φv + Ξ(ρ, ϕ, ψ, u, v, s, n)
    end
    function Ξ(
            ρ, ϕ, ψ, u, v, k, s, n, 
            n_a::SVector{N,T}, 
            t::SVector{N,T}) where {N, T<:Real}
        Φv = @. t * Φ(s*ρ, ϕ, n_a, α)
        @fastmath Ξv = k * sum(Φv) + Ξ(ρ, ϕ, ψ, u, v, s, n)
    end
    function ξx(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
        ret = ξx(ρ, ϕ, Ξv, s)
        return ret
    end
    function ξy(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
        ret = ξy(ρ, ϕ, Ξv, s)
        return ret
    end
    function ξz(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
        Ξv = Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)
        ret = ξz(ρ, ϕ, Ξv, s)
        return ret
    end
    function psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol = 1e-3, atol = 1e-4)
        resx = hcubature(
            dr -> ξx(dr[1], dr[2], ψ, u, v, k, s, α, n, n_a, t),
            [0.0,0.0],
            [1.0,2π],
            rtol = rtol,atol = atol)
        resy = hcubature(
            dr -> ξy(dr[1], dr[2], ψ, u, v, k, s, α, n, n_a, t),
                [0.0,0.0],
                [1.0,2π],
            rtol = rtol,atol = atol)
        resz = hcubature(
            dr -> ξz(dr[1], dr[2], ψ, u, v, k, s, α, n, n_a, t),
                [0.0,0.0],
                [1.0,2π],
            rtol = rtol,atol = atol)
        ret = 0.0
        ret += abs2(resx[1]::ComplexF64)
        ret += abs2(resy[1]::ComplexF64)
        ret += abs2(resz[1]::ComplexF64)
        return  ret/ pi^2
    end
    function PSF(
        x,y,z,
        λ,
        NA,
        plate::Nothing;
        rtol = 1e-3,
        atol = 1e-4)
        n = 1.
        s = NA #EP/R
        k = 2π / λ
        u = k * (1-sqrt(1- s^2))  / (1-sqrt(1- s^2/n^2)) * z
        v = k * s * sqrt(x^2 + y^2)
        ψ = atan(x,y)
        return psf_inner(ψ, u, v, s, n; rtol = rtol,atol = atol)
    end
    function PSF(
        x,y,z,
        λ,
        obj::Objective,
        plate::Nothing;
        rtol = 1e-3,
        atol = 1e-4)
        return PSF(
            x,y,z,
            λ,
            obj.NA,
            plate::Nothing;
            rtol = rtol,
            atol = atol)
    end
    function PSF(
            x,y,z,
            λ,
            NA,
            plate::PlaneParallelPlate;
            object::PlaneParallelPlate = plate,
            rtol = 1e-3,atol = 1e-4
            )
        n = object.n_λ(λ)
        n_a = plate.n_λ(λ)
        t = plate.t
        s = NA #EP/R
        k = 2π / λ
        u = k * (1-sqrt(1- s^2))  / (1-sqrt(1- s^2/n^2)) * z
        v = k * s * sqrt(x^2 + y^2)
        ψ = atan(x,y)
        return psf_inner(ψ, u, v, k, s, n, n_a, t; rtol = rtol,atol = atol)
    end
    function PSF(
        x,y,z,
        λ,
        obj::Objective,
        plate::PlaneParallelPlate;
        object::PlaneParallelPlate = plate,
        rtol = 1e-3,atol = 1e-4
        )
        return PSF(x,y,z,λ,obj.NA,plate; object = object, rtol = rtol, atol = atol)
    end
    function PSF(
            x,y,z,
            λ,
            NA,
            α,
            plate::PlaneParallelPlate;
            object::PlaneParallelPlate = plate,
            rtol = 1e-3,atol = 1e-4)
        n = object.n_λ(λ)
        n_a = plate.n_λ(λ)
        t = plate.t
        s = NA #EP/R
        k = 2π / λ
        u = k * (1-sqrt(1- s^2))  / (1-sqrt(1- s^2/n^2)) * z
        v = k * s * sqrt(x^2 + y^2)
        ψ = atan(x,y)
        return psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol = rtol,atol = atol)
    end
    function PSF(
        x,y,z,
        λ,
        obj::Objective,
        α,
        plate::PlaneParallelPlate;
        object::PlaneParallelPlate = plate,
        rtol = 1e-3,atol = 1e-4)
        n = object.n_λ(λ)
        n_a = plate.n_λ(λ)
        t = plate.t
        s = obj.NA #EP/R
        k = 2π / λ
        u = k * (1-sqrt(1- s^2))  / (1-sqrt(1- s^2/n^2)) * z
        v = k * s * sqrt(x^2 + y^2)
        ψ = atan(x,y)
        return psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol = rtol, atol = atol)
    end

    function PSF(
        x,y,z,
        λ,
        NA::Float64,
        plate::Nothing;
        object = plate,
        rtol = 1e-3,atol = 1e-4)
        if object == nothing
            n = 1.0
        else
            n = object.n_λ(λ)
        end
        s = NA #EP/R
        k = 2π / λ
        u = k * (1-sqrt(1- s^2))  / (1-sqrt(1- s^2/n^2)) * z
        v = k * s * sqrt(x^2 + y^2)
        ψ = atan(x,y)
        return psf_inner(ψ, u, v, s, n; rtol = rtol, atol = atol)
    end
    spl = SmoothingSplines.fit(SmoothingSplines.SmoothingSpline,wavelength,spectrum,250.0)
    function PSF_NV(x,y,z,λs::AbstractVector,obj::Objective, plate::Nothing;
        rtol = 1e-3,atol = 1e-4)
        NA = obj.NA
        w = SmoothingSplines.predict(spl,Vector(λs) .* 1000)
        weight = w/sum(w)
        res = PSF.(x,y,z,λs,NA,plate; rtol = rtol,atol = atol)
        return dot(res, weight)
    end
    function PSF_NV(x,y,z,λs::AbstractVector,obj::Objective, plate::PlaneParallelPlate; rtol = 1e-3,atol = 1e-4)
        NA = obj.NA
        w = SmoothingSplines.predict(spl,Vector(λs) .* 1000)
        weight = w/sum(w)
        res = map(λ -> PSF(x,y,z,λ,NA,plate;rtol = rtol,atol = atol), λs)
        return dot(res, weight)
    end
    function PSF_NV(x,y,z,λs::AbstractVector,obj::Objective, plate::PlaneParallelPlate, ro; rtol = 1e-3,atol = 1e-4)
        n = plate.n_λ(mean(λs))
        f = obj.f + (n-1) * t /n ## Paraxial Focal Shift
        NA = obj.NA
        α = -asin(ro / sqrt((f + z)^2 + x^2 + y^2))
        w = SmoothingSplines.predict(spl,Vector(λs) .* 1000)
        weight = w/sum(w)
        res = map(λ -> PSF(x,y,z,λ,NA,α, plate;rtol = rtol,atol = atol), λs)
        return dot(res, weight)
    end
    function Strehl(t,λ,obj::Objective)
        I_0 = PSF(0.,0.,0.,λ,obj,nothing)
        objective(z,t) = PSF(
            0,0,z,λ,obj,Diamond(t);rtol = 1e-10, atol = 1e-10) ./I_0
        res = optimize(z -> -objective(z,t), 0.0,2.0)
        return -minimum(res)
    end
    function Strehl(t,λ,α,obj::Objective)
        I_0 = PSF(0.,0.,0.,λ,obj,nothing)
        objective(z,t) = PSF(
            0,0,z,λ,
            obj,α,Diamond(t);rtol = 1e-10, atol = 1e-10) ./I_0
        res = optimize(z -> -objective(z,t), 0.0,2.0)
        return -minimum(res)
    end
    export PSF, PSF_NV, MPlanApo100x, Diamond
end
