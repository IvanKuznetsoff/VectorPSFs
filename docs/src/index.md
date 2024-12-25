# Overview

**VectorPSFs.jl** provides a comprehensive toolkit for **computing fully vectorial point spread functions (PSFs)** under various conditions, including the presence or absence of **aberrations** induced by plane-parallel plates. The package features specialized optimizations for modeling the PSFs of **NV centers** and supports **quantitative aberration analysis** by computing Strehl ratios and identifying plate thicknesses that preserve near-diffraction-limited performance ($\mathrm{Strehl} > 0.8$).

## Features

- **PSF Computations**:  
  The `PSF(...)` function handles scenarios involving normal or tilted incidence, immersion objectives, and spectral weighting. Users can specify `(x, y, z)` coordinates in micrometers, a chosen objective, and optionally a plane-parallel plate with a tilt angle.

- **Aberrations**:  
  Refractive-index mismatched transparent materials situated in the path of converging rays introduce **aberrations**. This package enables the computation of aberrated PSFs using a fully vectorial approach.

- **Plane-Parallel Plate Materials**:  
  Create `PlaneParallelPlate` objects for various materials, including diamond, fused silica, BK7 (borosilicate crown), magnesium fluoride, sapphire, or custom transparent plates whose refractive indices are given by Sellmeier equation.

- **Aberration Quantification**:  
  Compute **Strehl ratios** across defocus ranges or tilt scenarios, and automate the search for plate thicknesses that ensure near-diffraction-limited performance (\(\mathrm{Strehl} \geq 0.8\)).

- **Microscope Objective Conditioning**:  
  Define `Objective` structures to capture parameters such as focal length, numerical aperture, refractive index, and immersion medium (e.g., air, oil).

- **NV-Center Integration**:  
  Includes a specialized workflow for NV photoluminescence, utilizing spline-interpolated spectral data to compute polychromatic PSFs at different wavelengths.

## Basic Steps

1. **Define Microscope Objectives**:  
   Create `Objective` objects by specifying numerical aperture, immersion medium, and focal length.

2. **Construct Plane-Parallel Plates**:  
   Build `PlaneParallelPlate` instances (e.g., diamond, fused silica, BK7) using either built-in materials or custom Sellmeier models.

3. **Compute PSFs**:  
   Use `PSF(...)` to compute PSFs by specifying `(x, y, z)` in micrometers along with the relevant optical parameters (e.g., objective, plate). For compile-time dispatch, utilize the `PSFParams` parametric struct.

For installation and setup instructions, see [Installation](#Installation). To explore practical code examples, refer to [Examples](#Examples).

---


# Installation

To install **VectorPSFs.jl**, open the Julia package manager (Pkg REPL mode) and run:

```julia
] add https://github.com/IvanLuznetsoff/VectorPSFs.jl
```

The package requires:
- **HCubature.jl** for multi-dimensional integration,
- **StaticArrays.jl** for efficient small-array handling,
- **SmoothingSplines.jl** for NV-spectrum fitting,
- **Optim.jl** for Strehl ratio optimizations,
- plus standard libraries like **LinearAlgebra**, **Statistics**, etc.

Once installed, you can `using VectorPSFs` in your Julia session:
```julia
using VectorPSFs
```

Check your environment’s compatibility with **Julia ≥ 1.6** (recommended 1.10 or above).

---

# Examples

## PSF Calculation

The function `PSF(x, y, z, ...)` computes the on-axis or off-axis field intensity, integrating vectorially over the back pupil.  
You can choose among multiple overloads:

- **No plate, direct calculation** (simple, air):
  ```julia
  PSF(x, y, z, λ, NA; rtol=1e-3, atol=1e-4)
  ```
- **No plate, with `Objective`**:
  ```julia
  PSF(x, y, z, λ, obj::Objective; rtol=1e-3, atol=1e-4)
  ```
- **Normal-incidence plate, with `Objective`**:
  ```julia
  PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; ...)
  ```
- **Tilted-incidence plate, with `Objective`**:
  ```julia
  PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate, α; ...)
  ```

Here, `(x, y, z)` are in micrometers (µm).  
`λ` is the wavelength in µm, `NA` is numerical aperture, `α` is tilt angle in radians.

**Example** (normal incidence, no aberration):
```julia
using VectorPSFs

obj = MPlanApo50x()         # e.g. 50x objective
λ   = 0.632                 # 632 nm = 0.632 µm
```

To obtain a PSF map on a plane:
```julia
xs = range(-1,1,201)
ys = range(-1,1,201)
psf_map_xy = PSF.(x, y', 0.0, λ, obj)
```

## Strehl Ratio

The **Strehl ratio** is a measure of aberration severity, defined by the ratio of maximum intensity of the aberrated PSF to that of an unaberrated (diffraction-limited) PSF.  
**Functions**:
```julia
Strehl(t::Float64, λ::Float64, obj::Objective; ...)
Strehl(t::Float64, λ::Float64, α::Float64, obj::Objective; ...)
```
These calls do an internal defocus optimization.  
If `S > 0.8`, we often call it “diffraction-limited.”

**Example**:
```julia
using VectorPSFs, Optim

obj = MPlanApo100x()
λ   = 0.7         # 700 nm
s = Strehl(50.0, λ, obj, Diamond; zrange=[0.0, 5.0])  # Diamond thickness=50 µm
println("Strehl ratio = ", s)
```

## NV-Center Weighted PSF

We also offer a polychromatic weighted summation for NV emission:  
```julia
PSF(x, y, z, obj::Objective, nv::NVCenter; ...)
```
where `NVCenter` is a struct storing smoothed spectral weights. This integrates over multiple wavelengths.  

**Example**:
```julia
using VectorPSFs

nv = NVCenter(NVspectrum.wavelength)  # smoothing-based weighting
obj = MPlanApo100x()
plate_diamond = Diamond(100.0)
val_NV = PSF(0.0, 0.0, 3.0, obj, nv, plate_diamond)
println("NV-weighted PSF intensity = ", val_NV)
```

---

## Additional Topics

- **PlaneParallelPlate**: 
  Constructors like `Diamond(t)`, `FusedSilica(t)`, `BorosilicateCrown(t)`, etc., each define a wavelength-dependent refractive index using Sellmeier equations.
- **Objective**:
  Typical usage: `MPlanApo50x()`, `MPlanApo100x()`, etc., or define your own `Objective(f, NA, n)`.
- **ParametricPSF**:
  A more advanced method is to use parametric structs (e.g., `PSFParams{Mode}`) for compile-time dispatch of different incidence modes (no plate, normal, tilted).

For more details, consult the source code in `psf_core.jl`, `plane_parallel_plates.jl`, `objectives.jl`, etc.


Here’s the proofread version of your citation text:

---

# Citation

This package was developed as part of academic research. If you use it in your research, we would appreciate if you could cite our work ([arXiv:2402.14422](https://arxiv.org/abs/2402.14422)).