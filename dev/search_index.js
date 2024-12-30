var documenterSearchIndex = {"docs":
[{"location":"apiref/#API-Reference","page":"References","title":"API Reference","text":"","category":"section"},{"location":"apiref/","page":"References","title":"References","text":"Below is a reference for VectorPSFs.jl.   To see the source code or additional details, visit the corresponding file references (e.g., psf_core.jl, plane_parallel_plates.jl, objectives.jl, etc.).","category":"page"},{"location":"apiref/","page":"References","title":"References","text":"","category":"page"},{"location":"apiref/#Module-Index","page":"References","title":"Module Index","text":"","category":"section"},{"location":"apiref/","page":"References","title":"References","text":"PSF\nstrehl\nmaxtol_thick\nObjective\nPlaneParallelPlate\nPSFParams","category":"page"},{"location":"apiref/","page":"References","title":"References","text":"Modules = [VectorPSFs]\nPrivate = true","category":"page"},{"location":"apiref/#VectorPSFs.VectorPSFs","page":"References","title":"VectorPSFs.VectorPSFs","text":"module VectorPSFs\n\nThe primary module for computing vector PSFs (point spread functions) with or without aberrations caused by plane-parallel plates. Also supports NV-center weighted PSFs, Strehl ratio calculations, and lens objectives.\n\n\n\n\n\n","category":"module"},{"location":"apiref/#VectorPSFs.NVCenter","page":"References","title":"VectorPSFs.NVCenter","text":"NVCenter\n\nA struct for storing an NV emission spectrum, specifying:\n\nλs::Vector{Float64}: the sampling wavelengths in micrometers (µm)\nweight::Vector{Float64}: weighting (normalized spectral intensity)\n\nExample usage:\n\n# Construct from an existing array of wavelengths & weights\ncenter = NVCenter([0.632, 0.64, 0.65], [0.3, 0.5, 0.2])\n\n# Or use the fitted NV default (pulling from `nv_spectrum_spline`):\ncenter2 = NVCenter(0.5:0.01:0.7)\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.NVCenter-Tuple{Real, Real}","page":"References","title":"VectorPSFs.NVCenter","text":"NVCenter(λs::Real, weight::Real)\n\nConstructor that takes arrays of wavelengths λs and weights weight (must be same length). Converts them to Float64 and returns a new NVCenter.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.NVCenter-Tuple{Vector{Float64}}","page":"References","title":"VectorPSFs.NVCenter","text":"NVCenter(λs::Vector{Float64})\n\nConstructor that uses the global NV spectrum spline (nv_spectrum_spline)  to obtain weights for the given array of wavelengths λs. Weights are normalized (so their sum is 1).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.NoAberration","page":"References","title":"VectorPSFs.NoAberration","text":"struct NoAberration <: PSFMode end\n\nNo plate, no tilt scenario.\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.NormalIncidence","page":"References","title":"VectorPSFs.NormalIncidence","text":"struct NormalIncidence <: PSFMode end\n\nNormal incidence with a plate scenario.\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.Objective","page":"References","title":"VectorPSFs.Objective","text":"Objective\n\nA microscope objective characterized by:\n\nf::Float64: focal length (µm)\nNA::Float64: numerical aperture\nn::Float64: refractive index of immersion medium (e.g., 1.0 for air)\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.PSFMode","page":"References","title":"VectorPSFs.PSFMode","text":"abstract type PSFMode\n\nAbstract type for dispatching different PSF modes  (e.g. NoAberration, NormalIncidence, TiltedIncidence).\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.PSFParams","page":"References","title":"VectorPSFs.PSFParams","text":"mutable struct PSFParams{M<:PSFMode}\n\nHolds parameters for a certain PSF mode M.\n\nλ::Float64: wavelength in µm\nobj::Objective: the objective\nplate: either a PlaneParallelPlate or nothing\nα::Float64: tilt angle (radians)\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.PSFParams-Tuple{Float64, Objective, PlaneParallelPlate, Float64}","page":"References","title":"VectorPSFs.PSFParams","text":"PSFParams(λ::Float64, obj::Objective, plate::PlaneParallelPlate, α::Float64)\n\nConstruct a PSFParams{TiltedIncidence}, i.e. tilted incidence with plate. α in radians.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSFParams-Tuple{Float64, Objective, PlaneParallelPlate}","page":"References","title":"VectorPSFs.PSFParams","text":"PSFParams(λ::Float64, obj::Objective, plate::PlaneParallelPlate)\n\nConstruct a PSFParams{NormalIncidence}, i.e. normal incidence with plate.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSFParams-Tuple{Float64, Objective}","page":"References","title":"VectorPSFs.PSFParams","text":"PSFParams(λ::Float64, obj::Objective)\n\nConstruct a PSFParams{NoAberration}, i.e. no plate, no tilt.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PlaneParallelPlate","page":"References","title":"VectorPSFs.PlaneParallelPlate","text":"PlaneParallelPlate\n\nA PlaneParallelPlate has thickness t (in micrometers) and a refractive index function n_λ(λ), where λ is in micrometers. The field n_λ returns the refractive index at that wavelength.\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.TiltedIncidence","page":"References","title":"VectorPSFs.TiltedIncidence","text":"struct TiltedIncidence <: PSFMode end\n\nTilted incidence with a plate scenario.\n\n\n\n\n\n","category":"type"},{"location":"apiref/#VectorPSFs.BorosilicateCrown-Tuple{Any}","page":"References","title":"VectorPSFs.BorosilicateCrown","text":"BorosilicateCrown(t)\n\nCreates a PlaneParallelPlate of thickness t (µm) for BK7 (borosilicate crown glass). Uses the Schott Sellmeier coefficients:\n\nB = [1.040, 0.230, 1.010]\nC = [0.006, 0.020, 103.560]\n\nSource: refractiveindex.info BK7 page.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.CustomPlate-Tuple{Float64, Vector{Float64}, Vector{Float64}}","page":"References","title":"VectorPSFs.CustomPlate","text":"CustomPlate(t, B, C)\n\nCreate a PlaneParallelPlate of thickness t (µm), with a refractive index  defined by the Sellmeier coefficients B and C, where\n\nn^2 = 1 + ∑( Bᵢ * λ² / (λ² - Cᵢ) )\n\n(λ in micrometers).  B and C must be the same length.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.Diamond-Tuple{Any}","page":"References","title":"VectorPSFs.Diamond","text":"Diamond(t)\n\nApproximate dispersion for diamond with thickness t (in µm). λ is assumed in micrometers (µm), but the formula internally converts λ to nm  by multiplying by 1e3 for the Sellmeier-like expression.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.FusedSilica-Tuple{Any}","page":"References","title":"VectorPSFs.FusedSilica","text":"FusedSilica(t)\n\nCreates a PlaneParallelPlate of thickness t (µm) for Fused Silica (SiO₂), using the Malitson Sellmeier coefficients:\n\nB = [0.6961663, 0.4079426, 0.8974794]\nC = [0.0684043^2, 0.1162414^2, 9.896161^2]\n\nSource: Malitson data.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.LMPLFLN100XBD-Tuple{}","page":"References","title":"VectorPSFs.LMPLFLN100XBD","text":"LMPLFLN100XBD()\n\nReturns an Objective resembling the Olympus LMPLFLN 100x BD:\n\nf=1800 µm\nNA=0.8\nn=1.0 (air)\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.M10x-Tuple{}","page":"References","title":"VectorPSFs.M10x","text":"M10x()\n\nReturns an Objective resembling the Newport M10x:\n\nf=16500 µm\nNA=0.25\nn=1.0 (air)\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.MPlanApo100x-Tuple{}","page":"References","title":"VectorPSFs.MPlanApo100x","text":"MPlanApo100x()\n\nReturns an Objective resembling the Mitsutoyo MPlanApo 100x:\n\nf=2000 µm\nNA=0.7\nn=1.0 (air)\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.MPlanApo50x-Tuple{}","page":"References","title":"VectorPSFs.MPlanApo50x","text":"MPlanApo50x()\n\nReturns an Objective resembling the Mitsutoyo MPlanApo 50x:\n\nf=4000 µm\nNA=0.55\nn=1.0 (air)\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.MagnesiumFluoride-Tuple{Any}","page":"References","title":"VectorPSFs.MagnesiumFluoride","text":"MagnesiumFluoride(t)\n\nCreates a PlaneParallelPlate of thickness t (µm) for MgF₂  using the ordinary-ray Sellmeier coefficients (Li):\n\nB = [0.48755108, 0.39875031, 2.3120353]\nC = [0.04338408^2, 0.09461442^2, 23.793604^2]\n\nSource: MgF₂ (Li) data.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-NTuple{5, Any}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, λ, NA; rtol=1e-3, atol=1e-4)\n\nAberration-free PSF with no objective structure, using air (n=1).  \n\nλ in µm\nNA numeric aperture\n(x, y, z) in µm\nrtol, atol tolerances\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Any, Float64, PlaneParallelPlate, Any}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, λ, NA::Float64, plate::PlaneParallelPlate, α; ...)\n\nTilted incidence with no explicit Objective.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Any, Float64, PlaneParallelPlate}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, λ, NA::Float64, plate::PlaneParallelPlate; ...)\n\nNormal incidence with a plane-parallel plate, no explicit Objective.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Any, Objective, PlaneParallelPlate, Any}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate, α; ...)\n\nTilted incidence with plane + objective.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Any, Objective, PlaneParallelPlate}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; ...)\n\nNormal incidence with a plane-parallel plate and an explicit Objective.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Any, Objective}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, λ, obj::Objective; rtol=1e-3, atol=1e-4)\n\nAberration-free PSF with an Objective specifying immersion index, NA, etc.   No plate, normal incidence in air or immersion.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Objective, NVCenter, PlaneParallelPlate, Float64}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate, α::Float64; rtol=1e-3, atol=1e-4)\n\nCompute a polychromatic PSF for NV emission with a tilted-incidence plate.   Again, integrates over nv.λs and sums with nv.weight.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Objective, NVCenter, PlaneParallelPlate}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate; rtol=1e-3, atol=1e-4)\n\nCompute a polychromatic PSF for NV emission with a normal-incidence plate.   Broadcasts over nv.λs, calling single-wavelength PSF(...) for each λ, then weights with nv.weight.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, Objective, NVCenter}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, obj::Objective, nv::NVCenter; rtol=1e-3, atol=1e-4)\n\nCompute a polychromatic PSF for an NV center (no plate scenario).   Integrates over all wavelengths in nv.λs, calling the single-wavelength PSF(...) for each λ and weighting by nv.weight.\n\n(x, y, z) in µm\nobj is an Objective\nnv::NVCenter is the NV spectrum data\nrtol, atol are integration tolerances\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, PSFParams{VectorPSFs.NoAberration}}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, param::PSFParams{NoAberration}; rtol=1e-3, atol=1e-4)\n\nDispatch for the no-aberration case (NoAberration).  Calls PSF(x, y, z, λ, obj; rtol, atol).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, PSFParams{VectorPSFs.NormalIncidence}}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, param::PSFParams{NormalIncidence}; rtol=1e-3, atol=1e-4)\n\nDispatch for normal incidence with a plate.  Calls PSF(x, y, z, λ, obj, plate; rtol, atol).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.PSF-Tuple{Any, Any, Any, PSFParams{VectorPSFs.TiltedIncidence}}","page":"References","title":"VectorPSFs.PSF","text":"PSF(x, y, z, param::PSFParams{TiltedIncidence}; rtol=1e-3, atol=1e-4)\n\nDispatch for tilted incidence with a plate. Calls PSF(x, y, z, λ, obj, plate, α; rtol, atol).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.Sapphire-Tuple{Any}","page":"References","title":"VectorPSFs.Sapphire","text":"Sapphire(t)\n\nCreates a PlaneParallelPlate of thickness t (µm) for sapphire (Al₂O₃), using the ordinary-ray Sellmeier coefficients (Li):\n\nB = [1.4313493, 0.65054713, 5.3414021]\nC = [0.0726631^2, 0.1193242^2, 18.028251^2]\n\nSource: Sapphire (Li) data.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.UPLXAPO100X-Tuple{}","page":"References","title":"VectorPSFs.UPLXAPO100X","text":"UPLXAPO100X()\n\nReturns an Objective resembling the Olympus UPLXAPO100X:\n\nf=1800 µm\nNA=1.45\nn=1.515 (oil immersion)\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.maxtol_thick","page":"References","title":"VectorPSFs.maxtol_thick","text":"maxtol_thick(λ::Float64, obj::Objective; plate=Diamond, t_range=(0.0, 500.0), zrange=[-2.0,4.0], ...)\n\nFind thickness t in t_range that yields a Strehl ratio near 0.8.   No tilt is assumed.   We define:\n\ncost(t) = (strehl(t, λ, obj, plate; zrange=zrange) - 0.8)^2\n\nand minimize with a 1D search.   Returns the thickness that best satisfies S ≈ 0.8.\n\n\n\n\n\n","category":"function"},{"location":"apiref/#VectorPSFs.maxtol_thick-Tuple{Float64, Objective, Function, Float64}","page":"References","title":"VectorPSFs.maxtol_thick","text":"maxtol_thick(λ::Float64, obj::Objective, plate::Function, α::Float64; t_range=(0.0, 500.0), zrange=[-2.0,4.0], ...)\n\nSame as above, but for tilted incidence angle α.  \n\ncost(t) = (strehl(t, λ, obj, plate, α; zrange=zrange) - 0.8)^2\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.psf_inner-NTuple{5, Any}","page":"References","title":"VectorPSFs.psf_inner","text":"psf_inner(ψ, u, v, s, n; rtol=1e-3, atol=1e-4)\n\nComputes the aberration-free PSF intensity by numerically integrating the x-, y-, and z-polarized integrands via hcubature.\n\nψ, u, v are phase parameters relating to off-axis, defocus, etc.\ns is the normalized aperture.\nn is the refractive index (1.0 or immersion).\nrtol, atol are integration tolerances.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.psf_inner-NTuple{8, Any}","page":"References","title":"VectorPSFs.psf_inner","text":"psf_inner(ψ, u, v, k, s, n, n_a, t; rtol=1e-3, atol=1e-4)\n\nCompute the PSF with normal incidence on a plane-parallel plate via hcubature.\n\nk is the wave number (2πn_obj / λ).\nt is plate thickness in µm.\nn_a is the plate’s refractive index (possibly scaled by immersion factor).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.psf_inner-NTuple{9, Any}","page":"References","title":"VectorPSFs.psf_inner","text":"psf_inner(ψ, u, v, k, s, α, n, n_a, t; rtol=1e-3, atol=1e-4)\n\nCompute the PSF with tilted incidence on a plane-parallel plate via hcubature.\n\nα is the tilt angle (radians).\nk is wave number (2π n_imm / λ).\nn_a is plate index scaled if needed.\nOther parameters as before.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.sellmeier-Tuple{Any, Vector{Float64}, Vector{Float64}}","page":"References","title":"VectorPSFs.sellmeier","text":"sellmeier(λ, B, C)\n\nGeneric Sellmeier equation:\n\nn^2 = 1 + ∑( Bᵢ * λ² / (λ² - Cᵢ) )\n\nλ in micrometers (µm).\nB/C are matching-length vectors of Sellmeier coefficients.\n\nReturns the refractive index n(λ).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.strehl","page":"References","title":"VectorPSFs.strehl","text":"strehl(t, λ, obj::Objective, plate::Function=Diamond; \n       zrange=[-2.0, 4.0], rtol=1e-10, atol=1e-10)\n\nCompute the Strehl ratio when a plate of thickness t (µm) is placed in the optical path.\n\nλ is wavelength in µm\nobj is the objective (e.g. MPlanApo50x())\nplate is a constructor function for a PlaneParallelPlate (defaults to Diamond).\nzrange in µm is the defocus search range: the function does a 1D optimization to find the best focus.\nrtol, atol are integration tolerances for the PSF.\n\nReturns the maximum on-axis intensity (aberrated) normalized by the unaberrated case. Example:\n\nS = strehl(50.0, 0.7, MPlanApo100x())  # diamond thickness=50 µm, λ=0.7 µm\n\n\n\n\n\n","category":"function"},{"location":"apiref/#VectorPSFs.strehl-2","page":"References","title":"VectorPSFs.strehl","text":"strehl(t, λ, obj::Objective, z_est::Float64, plate::Function=Diamond; rtol=1e-10, atol=1e-10)\n\nCompute the on-axis Strehl ratio for a plate thickness t (µm) at a specified defocus z_est (µm). No tilt assumed. \n\n\n\n\n\n","category":"function"},{"location":"apiref/#VectorPSFs.strehl-Tuple{Float64, Float64, Objective, Float64, Function, Float64}","page":"References","title":"VectorPSFs.strehl","text":"strehl(t, λ, obj::Objective, z_est::Float64, plate::Function, α::Float64; rtol=1e-10, atol=1e-10)\n\nSame as above, but with incidence angle α (tilt).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.strehl-Tuple{Float64, Float64, Objective, Function, Float64}","page":"References","title":"VectorPSFs.strehl","text":"strehl(t, λ, obj::Objective, plate::Function, α::Float64; \n       zrange=[-2.0, 4.0], rtol=1e-10, atol=1e-10)\n\nCompute the Strehl ratio for a plate of thickness t (µm), wavelength λ (µm),  Objective obj, and tilt angle α (radians).   The search for best defocus is done in zrange (µm).  \n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.Ξ-NTuple{10, Any}","page":"References","title":"VectorPSFs.Ξ","text":"Ξ(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)\n\nPhase term for normal incidence with a plane-parallel plate. Adds the plate-induced phase k * t * Φ(...) to the aberration-free phase.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.Ξ-NTuple{11, Any}","page":"References","title":"VectorPSFs.Ξ","text":"Ξ(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)\n\nPhase term for tilted incidence with a plane-parallel plate. Adds the tilt-based phase Φ(..., α) plus the base aberration-free term.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.Ξ-NTuple{7, Any}","page":"References","title":"VectorPSFs.Ξ","text":"Ξ(ρ, ϕ, ψ, u, v, s, n)\n\nComputes the phase term for the aberration-free case.  \n\nρ, ϕ are polar coordinates in the pupil.\nψ is an azimuthal shift (e.g. relating to off-axis propagation).\nu, v are defocus/piston terms in the exponent.\ns is the normalized aperture (NA).\nn is the refractive index (often 1.0 for air or obj.n for immersion).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.Φ-NTuple{4, Any}","page":"References","title":"VectorPSFs.Φ","text":"Φ(sρ, ϕ, n, α)\n\nComputes the aberration phase for tilted incidence.\n\nsρ: again, s * ρ.\nϕ: azimuthal angle in the pupil plane.\nn: refractive index.\nα: tilt angle (radians).\n\nReturns a phase shift accounting for tilt plus the normal-incidence phase.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.Φ-Tuple{Any, Any}","page":"References","title":"VectorPSFs.Φ","text":"Φ(sρ, n)\n\nComputes the aberration phase for normal incidence.\n\nsρ: the product of normalized aperture s and radial coordinate ρ in the pupil plane.\nn: refractive index (e.g., 1.0 for air or the immersion index).\n\nReturns a dimensionless phase shift due to the optical path difference.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξx-NTuple{10, Any}","page":"References","title":"VectorPSFs.ξx","text":"ξx(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)\n\nx-polarized integrand with normal-incidence plate. \n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξx-NTuple{11, Any}","page":"References","title":"VectorPSFs.ξx","text":"ξx(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)\n\nx-polarized integrand with tilted-incidence plate.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξx-NTuple{4, Any}","page":"References","title":"VectorPSFs.ξx","text":"ξx(ρ, ϕ, Ξ, s)\n\nx-polarized integrand for the aberration-free PSF integral.  \n\nρ, ϕ are pupil-plane coordinates.\nΞ is the accumulated phase from Ξ(...).\ns is the normalized aperture.\n\nReturns the complex field contribution in the x-direction.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξy-NTuple{10, Any}","page":"References","title":"VectorPSFs.ξy","text":"ξy(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)\n\ny-polarized integrand with normal-incidence plate.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξy-NTuple{11, Any}","page":"References","title":"VectorPSFs.ξy","text":"ξy(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)\n\ny-polarized integrand with tilted-incidence plate.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξy-NTuple{4, Any}","page":"References","title":"VectorPSFs.ξy","text":"ξy(ρ, ϕ, Ξ, s)\n\ny-polarized integrand for the aberration-free PSF integral. Similar arguments as ξx(...).\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξz-NTuple{10, Any}","page":"References","title":"VectorPSFs.ξz","text":"ξz(ρ, ϕ, ψ, u, v, k, s, n, n_a, t)\n\nz-polarized integrand with normal-incidence plate.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξz-NTuple{11, Any}","page":"References","title":"VectorPSFs.ξz","text":"ξz(ρ, ϕ, ψ, u, v, k, s, α, n, n_a, t)\n\nz-polarized integrand with tilted-incidence plate.\n\n\n\n\n\n","category":"method"},{"location":"apiref/#VectorPSFs.ξz-NTuple{4, Any}","page":"References","title":"VectorPSFs.ξz","text":"ξz(ρ, ϕ, Ξ, s)\n\nz-polarized integrand for the aberration-free PSF integral. Similar arguments as ξx(...).\n\n\n\n\n\n","category":"method"},{"location":"apiref/","page":"References","title":"References","text":"Modules = [VectorPSFs.NVspectrum]\nPrivate = true","category":"page"},{"location":"apiref/#VectorPSFs.NVspectrum","page":"References","title":"VectorPSFs.NVspectrum","text":"module NVspectrum\n\nNV center spectrum. \n\nThe wavelength array is given in micrometers (µm).\nThe spectrum array holds corresponding intensities (arbitrary units).\nnv_spectrum_spline is a smoothing spline fit of spectrum vs. wavelength.\n\nYou can access the raw data via NVspectrum.wavelength and NVspectrum.spectrum, or use nv_spectrum_spline for continuous interpolation.\n\n\n\n\n\n","category":"module"},{"location":"apiref/#VectorPSFs.NVspectrum.nv_spectrum_spline","page":"References","title":"VectorPSFs.NVspectrum.nv_spectrum_spline","text":"nv_spectrum_spline\n\nA smoothing spline fit (SmoothingSplines.SmoothingSpline) constructed from  wavelength and spectrum, using a smoothing parameter of 250.0. This allows continuous interpolation/extrapolation of the NV emission data.\n\n\n\n\n\n","category":"constant"},{"location":"apiref/#VectorPSFs.NVspectrum.spectrum","page":"References","title":"VectorPSFs.NVspectrum.spectrum","text":"spectrum\n\nA constant array of raw intensities (arbitrary units) for each entry in wavelength.\n\n\n\n\n\n","category":"constant"},{"location":"apiref/#VectorPSFs.NVspectrum.wavelength","page":"References","title":"VectorPSFs.NVspectrum.wavelength","text":"wavelength\n\nA constant array of NV emission wavelengths (in µm).\n\n\n\n\n\n","category":"constant"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"To install VectorPSFs.jl, open the Julia package manager (Pkg REPL mode) and run:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"] add https://github.com/IvanLuznetsoff/VectorPSFs.jl","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"The package requires:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"HCubature.jl for multi-dimensional integration,\nStaticArrays.jl for efficient small-array handling,\nSmoothingSplines.jl for NV-spectrum fitting,\nOptim.jl for Strehl ratio optimizations,\nplus standard libraries like LinearAlgebra, Statistics, etc.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Once installed, you can using VectorPsfs in your Julia session:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"using VectorPsfs","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Check your environment’s compatibility with Julia ≥ 1.6 (recommended 1.10 or above). ```","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"VectorPSFs.jl provides a comprehensive toolkit for computing fully vectorial point spread functions (PSFs) under various conditions, including the presence or absence of aberrations induced by plane-parallel plates. The package features specialized optimizations for modeling the PSFs of NV centers and supports quantitative aberration analysis by computing Strehl ratios and identifying plate thicknesses that preserve near-diffraction-limited performance (mathrmStrehl  08).","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PSF Computations:   The PSF(...) function handles scenarios involving normal or tilted incidence, immersion objectives, and spectral weighting. Users can specify (x, y, z) coordinates in micrometers, a chosen objective, and optionally a plane-parallel plate with a tilt angle.\nAberrations:   Refractive-index mismatched transparent materials situated in the path of converging rays introduce aberrations. This package enables the computation of aberrated PSFs using a fully vectorial approach.\nPlane-Parallel Plate Materials:   Create PlaneParallelPlate objects for various materials, including diamond, fused silica, BK7 (borosilicate crown), magnesium fluoride, sapphire, or custom transparent plates whose refractive indices are given by Sellmeier equation.\nAberration Quantification:   Compute Strehl ratios across defocus ranges or tilt scenarios, and automate the search for plate thicknesses that ensure near-diffraction-limited performance ((\\mathrm{Strehl} \\geq 0.8)).\nMicroscope Objective Conditioning:   Define Objective structures to capture parameters such as focal length, numerical aperture, refractive index, and immersion medium (e.g., air, oil).\nNV-Center Integration:   Includes a specialized workflow for NV photoluminescence, utilizing spline-interpolated spectral data to compute polychromatic PSFs at different wavelengths.","category":"page"},{"location":"#Basic-Steps","page":"Home","title":"Basic Steps","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Define Microscope Objectives:   Create Objective objects by specifying numerical aperture, immersion medium, and focal length.\nConstruct Plane-Parallel Plates:   Build PlaneParallelPlate instances (e.g., diamond, fused silica, BK7) using either built-in materials or custom Sellmeier models.\nCompute PSFs:   Use PSF(...) to compute PSFs by specifying (x, y, z) in micrometers along with the relevant optical parameters (e.g., objective, plate). For compile-time dispatch, utilize the PSFParams parametric struct.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For installation and setup instructions, see Installation. To explore practical code examples, refer to Examples.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install VectorPSFs.jl, open the Julia package manager (Pkg REPL mode) and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add https://github.com/IvanLuznetsoff/VectorPSFs.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package requires:","category":"page"},{"location":"","page":"Home","title":"Home","text":"HCubature.jl for multi-dimensional integration,\nStaticArrays.jl for efficient small-array handling,\nSmoothingSplines.jl for NV-spectrum fitting,\nOptim.jl for Strehl ratio optimizations,\nplus standard libraries like LinearAlgebra, Statistics, etc.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Once installed, you can using VectorPSFs in your Julia session:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using VectorPSFs","category":"page"},{"location":"","page":"Home","title":"Home","text":"Check your environment’s compatibility with Julia ≥ 1.6 (recommended 1.10 or above).","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"#PSF-Calculation","page":"Home","title":"PSF Calculation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The function PSF(x, y, z, ...) computes the on-axis or off-axis field intensity, integrating vectorially over the back pupil.   You can choose among multiple overloads:","category":"page"},{"location":"","page":"Home","title":"Home","text":"No plate, direct calculation (simple, air):\nPSF(x, y, z, λ, NA; rtol=1e-3, atol=1e-4)\nNo plate, with Objective:\nPSF(x, y, z, λ, obj::Objective; rtol=1e-3, atol=1e-4)\nNormal-incidence plate, with Objective:\nPSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; ...)\nTilted-incidence plate, with Objective:\nPSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate, α; ...)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here, (x, y, z) are in micrometers (µm).   λ is the wavelength in µm, NA is numerical aperture, α is tilt angle in radians.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example (normal incidence, no aberration):","category":"page"},{"location":"","page":"Home","title":"Home","text":"using VectorPSFs\n\nobj = MPlanApo50x()         # e.g. 50x objective\nλ   = 0.632                 # 632 nm = 0.632 µm","category":"page"},{"location":"","page":"Home","title":"Home","text":"To obtain a PSF map on a plane:","category":"page"},{"location":"","page":"Home","title":"Home","text":"xs = range(-1,1,201)\nys = range(-1,1,201)\npsf_map_xy = PSF.(x, y', 0.0, λ, obj)","category":"page"},{"location":"#Strehl-Ratio","page":"Home","title":"Strehl Ratio","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Strehl ratio is a measure of aberration severity, defined by the ratio of maximum intensity of the aberrated PSF to that of an unaberrated (diffraction-limited) PSF.   Functions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Strehl(t::Float64, λ::Float64, obj::Objective; ...)\nStrehl(t::Float64, λ::Float64, α::Float64, obj::Objective; ...)","category":"page"},{"location":"","page":"Home","title":"Home","text":"These calls do an internal defocus optimization.   If S > 0.8, we often call it “diffraction-limited.”","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using VectorPSFs, Optim\n\nobj = MPlanApo100x()\nλ   = 0.7         # 700 nm\ns = Strehl(50.0, λ, obj, Diamond; zrange=[0.0, 5.0])  # Diamond thickness=50 µm\nprintln(\"Strehl ratio = \", s)","category":"page"},{"location":"#NV-Center-Weighted-PSF","page":"Home","title":"NV-Center Weighted PSF","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We also offer a polychromatic weighted summation for NV emission:  ","category":"page"},{"location":"","page":"Home","title":"Home","text":"PSF(x, y, z, obj::Objective, nv::NVCenter; ...)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where NVCenter is a struct storing smoothed spectral weights. This integrates over multiple wavelengths.  ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using VectorPSFs\n\nnv = NVCenter(NVspectrum.wavelength)  # smoothing-based weighting\nobj = MPlanApo100x()\nplate_diamond = Diamond(100.0)\nval_NV = PSF(0.0, 0.0, 3.0, obj, nv, plate_diamond)\nprintln(\"NV-weighted PSF intensity = \", val_NV)","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Additional-Topics","page":"Home","title":"Additional Topics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PlaneParallelPlate:  Constructors like Diamond(t), FusedSilica(t), BorosilicateCrown(t), etc., each define a wavelength-dependent refractive index using Sellmeier equations.\nObjective: Typical usage: MPlanApo50x(), MPlanApo100x(), etc., or define your own Objective(f, NA, n).\nParametricPSF: A more advanced method is to use parametric structs (e.g., PSFParams{Mode}) for compile-time dispatch of different incidence modes (no plate, normal, tilted).","category":"page"},{"location":"","page":"Home","title":"Home","text":"For more details, consult the source code in psf_core.jl, plane_parallel_plates.jl, objectives.jl, etc.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here’s the proofread version of your citation text:","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Citation","page":"Home","title":"Citation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package was developed as part of academic research. If you use it in your research, we would appreciate if you could cite our work (arXiv:2402.14422).","category":"page"},{"location":"example/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"example/#PSF-Calculation","page":"Examples","title":"PSF Calculation","text":"","category":"section"},{"location":"example/","page":"Examples","title":"Examples","text":"The function PSF(x, y, z, ...) computes the on-axis or off-axis field intensity, integrating vectorially over the back pupil.   You can choose among multiple overloads:","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"No plate, direct calculation (simple, air):\nPSF(x, y, z, λ, NA; rtol=1e-3, atol=1e-4)\nNo plate, with Objective:\nPSF(x, y, z, λ, obj::Objective; rtol=1e-3, atol=1e-4)\nNormal-incidence plate, with Objective:\nPSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; ...)\nTilted-incidence plate, with Objective:\nPSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate, α; ...)","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"Here, (x, y, z) are in micrometers (µm).   λ is the wavelength in µm, NA is numerical aperture, α is tilt angle in radians.","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"Example (normal incidence, no aberration):","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"using VectorPsfs\n\nobj = MPlanApo50x()         # e.g. 50x objective\nλ   = 0.632                 # 632 nm = 0.632 µm","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"To obtain a PSF map on a plane:","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"xs = range(-1,1,201)\nys = range(-1,1,201)\npsf_map_xy = [PSF(x, y, 0.0, λ, obj) for x in xs, y in ys]","category":"page"},{"location":"example/#Strehl-Ratio","page":"Examples","title":"Strehl Ratio","text":"","category":"section"},{"location":"example/","page":"Examples","title":"Examples","text":"The Strehl ratio is a measure of aberration severity, defined by the ratio of maximum intensity of the aberrated PSF to that of an unaberrated (diffraction-limited) PSF.   Functions:","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"Strehl(t::Float64, λ::Float64, obj::Objective; ...)\nStrehl(t::Float64, λ::Float64, α::Float64, obj::Objective; ...)","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"These calls do an internal defocus optimization.   If S > 0.8, we often call it “diffraction-limited.”","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"Example:","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"using VectorPsfs, Optim\n\nobj = MPlanApo100x()\nλ   = 0.7         # 700 nm\ns = Strehl(50.0, λ, obj, Diamond; zrange=[0.0, 5.0])  # Diamond thickness=50 µm\nprintln(\"Strehl ratio = \", s)","category":"page"},{"location":"example/#NV-Center-Weighted-PSF","page":"Examples","title":"NV-Center Weighted PSF","text":"","category":"section"},{"location":"example/","page":"Examples","title":"Examples","text":"We also offer a polychromatic approach for NV emission:  ","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"PSF(x, y, z, obj::Objective, nv::NVCenter; ...)","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"where NVCenter is a struct storing smoothed spectral weights. This integrates over multiple wavelengths.  ","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"Example:","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"using VectorPsfs\n\nnv = NVCenter(NVspectrum.wavelength)  # smoothing-based weighting\nobj = MPlanApo100x()\nplate_diamond = Diamond(100.0)\nval_NV = PSF(0.0, 0.0, 3.0, obj, nv, plate_diamond)\nprintln(\"NV-weighted PSF intensity = \", val_NV)","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"example/#Additional-Topics","page":"Examples","title":"Additional Topics","text":"","category":"section"},{"location":"example/","page":"Examples","title":"Examples","text":"PlaneParallelPlate:  Constructors like Diamond(t), FusedSilica(t), BorosilicateCrown(t), etc., each define a wavelength-dependent refractive index using Sellmeier equations.\nObjective: Typical usage: MPlanApo50x(), MPlanApo100x(), etc., or define your own Objective(f, NA, n).\nParametricPSF: A more advanced method is to use parametric structs (e.g., PSFParams{Mode}) for compile-time dispatch of different incidence modes (no plate, normal, tilted).","category":"page"},{"location":"example/","page":"Examples","title":"Examples","text":"For more details, consult the source code in psf_core.jl, plane_parallel_plates.jl, objectives.jl, etc. ```","category":"page"}]
}
