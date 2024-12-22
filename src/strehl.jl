# strehl.jl

import ..psf_core: PSF
import ..plane_parallel_plates: Diamond
import ..objectives: Objective
import ..nvcenter: NVCenter  # if needed

function Strehl(t, λ, obj::Objective)
    I_0 = PSF(0.0, 0.0, 0.0, λ, obj, nothing)
    objective_func(z) = PSF(0.0, 0.0, z, λ, obj, Diamond(t); rtol=1e-10, atol=1e-10) / I_0
    res = optimize(z -> -objective_func(z), 0.0, 2.0)
    return -minimum(res)
end

function Strehl(t, λ, α, obj::Objective)
    I_0 = PSF(0.0, 0.0, 0.0, λ, obj, nothing)
    objective_func(z) = PSF(0.0, 0.0, z, λ, obj, α, Diamond(t); rtol=1e-10, atol=1e-10) / I_0
    res = optimize(z -> -objective_func(z), 0.0, 2.0)
    return -minimum(res)
end
