# Home

## ToC

```@contents
Pages = [
    "installation.md"
    "example.md"
]
Depth = 1
```

---

## Overview

`VectorPSFs.jl` provides tools for computing **vectorial point spread functions (PSFs)**, with or without aberrations, caused by ray transmission through a plane-parallel plate. 
It is particularly optimized for calculating the PSFs of **NV centers**.
The package also enables the quantification of aberration effects by **computing Strehl ratios** and determining plate thicknesses that achieve the diffraction limit (Strehl ratio > 0.8). 

For PSF conditionings, follow these steps:
1. **Define microscope objectives** using `Objective` objects, specifying parameters such as numerical aperture, immersion medium, and focal length.
2. **Construct plane-parallel plates** (e.g., diamond) using Sellmeier dispersion models.


---

## Submodules and Key Components

### 1. Objectives

The file `objectives.jl` defines a simple `Objective` struct that encapsulates:

- `f`: focal length (e.g. 2000 μm),
- `NA`: numerical aperture (dimensionless),
- `n`: refractive index in the immersion medium (e.g. 1.0 for air, 1.515 for oil, etc.).

Some standard objectives are provided as convenience constructors:

- `MPlanApo100x()`
- `LMPLFLN100XBD()`
- `MPlanApo50x()`
- `M10x()`
- `UPLXAPO100X()`

**Example Usage**:
```julia
using VectorPSFs

obj = MPlanApo100x()  # f=2000, NA=0.7, n=1.0
```

### 2. Plane-Parallel Plates

The file `plane_parallel_plates.jl` provides:

- A `PlaneParallelPlate` struct:  
  - `t` (thickness in micrometers, μm)
  - `n_λ` (a function giving refractive index vs. wavelength λ in micrometers)

Predefined materials using Sellmeier equations include:

- `Diamond(t)`, `FusedSilica(t)`, `BorosilicateCrown(t)`, `Sapphire(t)`, `MagnesiumFluoride(t)`, and a `CustomPlate(t, B, C)` for user-supplied Sellmeier coefficients.

**Example Usage**:
```julia
plate_diamond = Diamond(100.0)  # 100 μm-thick diamond plate
plate_bk7     = BorosilicateCrown(250.0)  # BK7 plate, 250 μm thick
```

### 3. PSF Computations

#### 3.1 High-level `PSF` functions

`PSF(x, y, z, ...)` is defined in multiple overloads:

- **No plate, no objective** (simple, air):
  ```julia
  PSF(x, y, z, λ, NA; rtol=1e-3, atol=1e-4)
  ```
- **No plate, with objective** (aberration-free, immersion):
  ```julia
  PSF(x, y, z, λ, obj::Objective; rtol=1e-3, atol=1e-4)
  ```
- **Normal-incidence plate, no objective**:
  ```julia
  PSF(x, y, z, λ, NA, plate::PlaneParallelPlate; ...)
  ```
- **Normal-incidence plate, with objective**:
  ```julia
  PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; ...)
  ```
- **Tilted-incidence plate, no objective**:
  ```julia
  PSF(x, y, z, λ, NA, plate::PlaneParallelPlate, α; ...)
  ```
- **Tilted-incidence plate, with objective**:
  ```julia
  PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate, α; ...)
  ```

`(x, y, z)` are spatial coordinates (in μm), `λ` is wavelength in micrometers, `NA` is the numerical aperture, `α` is the tilt angle in radians.

The various `rtol`/`atol` arguments control the tolerance in the numerical integration (via `HCubature`).

**Example** (Normal incidence, plate, objective):
```julia
using VectorPSFs

obj = MPlanApo50x()           # e.g. 50x objective
plate = Diamond(50.0)         # 50 μm diamond
λ = 0.632                     # 632 nm = 0.632 μm
psf_val = PSF(0.0, 0.0, 0.0, λ, obj, plate)  
println("PSF on-axis intensity: ", psf_val)
```

#### 3.2 `PSFParams` Structures

`PSFParams` wraps the relevant parameters (wavelength, objective, plate, tilt) along with a type marker `PSFMode` indicating whether it is `NoAberration`, `NormalIncidence`, or `TiltedIncidence`. 

You can construct it like so:
- **No plate**:
  ```julia
  param = PSFParams(0.632, MPlanApo100x())
  I = PSF(x, y, z, param)
  ```
- **Normal incidence**:
  ```julia
  param = PSFParams(0.532, MPlanApo50x(), Diamond(100.0))
  I = PSF(0.0, 0.0, 1.0, param)
  ```
- **Tilted incidence**:
  ```julia
  param = PSFParams(0.532, MPlanApo50x(), Diamond(100.0), α=0.1)
  I = PSF(0.0, 0.0, 1.0, param)
  ```

### 4. NV-Center Emission Spectrum

In `NVspectrum.jl`, the module defines arrays of `wavelength` and corresponding raw `spectrum` data, then creates a smoothing spline (`nv_spectrum_spline`).

An `NVCenter` struct stores:
```julia
struct NVCenter
    λs::Vector{Float64}   
    weight::Vector{Float64}
end
```
and you can construct it in multiple ways. A typical usage is:
```julia
using VectorPSFs

nv_default = NVCenter(NVspectrum.wavelength)  # uses the default spline
```
Then, you can compute the polychromatic PSF from an NV emitter using:
```julia
PSF(x, y, z, obj::Objective, nv::NVCenter; rtol, atol)
PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate; ...)
PSF(x, y, z, obj::Objective, nv::NVCenter, plate::PlaneParallelPlate, α; ...)
```
It integrates over all wavelengths in `nv.λs`, weighting each by `nv.weight`.

**Example** (NV center with objective and diamond plate):
```julia
using VectorPSFs

obj       = MPlanApo100x()
plate     = Diamond(100.0)
nvcenter  = NVCenter(NVspectrum.wavelength)  # or simply NVCenter(NVspectrum.wavelength)
x, y, z   = 0.0, 0.0, 0.0
psf_intensity = PSF(x, y, z, obj, nvcenter, plate)
println("Weighted PSF from NV-center emission = ", psf_intensity)
```

### 5. Strehl Ratio Computations

The `strehl.jl` file provides `strehl(...)` methods that:

1. Compute on-axis intensity (x=0, y=0) as a function of defocus `z`.
2. Find the maximum intensity in a given `zrange`.
3. Normalize by the unaberrated (or reference) case to get a Strehl ratio.

**Method signatures**:

```julia
strehl(t::Float64, λ::Float64, obj::Objective;
       plate::Function=Diamond, zrange=[-2.0, 4.0], rtol=1e-10, atol=1e-10)

strehl(t::Float64, λ::Float64, obj::Objective, plate::Function, α::Float64;
       zrange=[-2.0, 4.0], ...)

strehl(t::Float64, λ::Float64, obj::Objective, z_est::Float64;
       plate::Function=Diamond, ...)

strehl(t::Float64, λ::Float64, obj::Objective, z_est::Float64,
       plate::Function, α::Float64; ...)
```

Here,
- `t` is plate thickness,
- `λ` is wavelength (μm),
- `obj` is the `Objective`,
- `α` is tilt angle (optional),
- `z_est` is a user-chosen defocus for direct evaluation,
- `zrange` is a bracket around which we do a 1D optimization for the best focus.

**Example**:
```julia
using VectorPSFs, Optim

obj = UPLXAPO100X()         # e.g. NA=1.45, n=1.515
λ   = 0.637                 # 637 nm
# Find the Strehl ratio for a 50 μm diamond plate:
S = strehl(50.0, λ, obj, Diamond; zrange=[-5.0, 5.0])
println("Strehl ratio = ", S)
```

### 6. Finding Plate Thickness for Desired Strehl

`maxtol_thick(...)` searches for a thickness `t` in a given range such that the Strehl ratio is as close to 0.8 (or any built-in target) as possible.

Two variants exist (normal incidence vs. tilted):
```julia
maxtol_thick(λ::Float64, obj::Objective; plate::Function=Diamond, t_range=(0.0, 500.), ...)
maxtol_thick(λ::Float64, obj::Objective, plate::Function, α::Float64; t_range=(0.0, 500.), ...)
```

**Example**:
```julia
using VectorPSFs, Optim

obj = MPlanApo100x()
λ   = 0.637
t_star = maxtol_thick(λ, obj, Diamond; t_range=(0.0, 200.0))
println("Thickness achieving Strehl ≈ 0.8: ", t_star, " μm")
```

---

## Example Workflows

Below is a more complete sequence showcasing key steps:

```julia
using VectorPSFs, Optim

# 1. Define your objective and plate
obj    = MPlanApo100x()          # 100x objective, NA=0.7, n=1.0
plate  = Diamond(50.0)           # 50 μm diamond

# 2. Choose wavelength (in μm)
λ_red  = 0.637                   # 637 nm

# 3. Compute on-axis PSF with/without plate
psf_no_plate     = PSF(0.0, 0.0, 0.0, λ_red, obj)           # no plate
psf_with_plate   = PSF(0.0, 0.0, 0.0, λ_red, obj, plate)    # normal-incidence diamond

println("On-axis, no plate:     ", psf_no_plate)
println("On-axis, with plate:   ", psf_with_plate)

# 4. Evaluate PSF off-axis
x, y, z = 0.5, 0.5, 1.0
val_off_axis = PSF(x, y, z, λ_red, obj, plate)
println("PSF at (x,y,z) = (0.5, 0.5, 1.0): ", val_off_axis)

# 5. Compute Strehl ratio for t=50 μm diamond
S_50 = strehl(50.0, λ_red, obj, Diamond; zrange=[-2.0, 2.0])
println("Strehl ratio at t=50 μm: ", S_50)

# 6. Find thickness that yields ~0.8 Strehl
t_star = maxtol_thick(λ_red, obj, Diamond; t_range=(0.0, 200.0))
println("Thickness achieving ~0.8 Strehl: ", t_star)
```

---

## NV-Center Example

A polychromatic PSF, summing over the NV emission:

```julia
using VectorPSFs

# Create an NVCenter (already has an internal spline from 460 nm to ~810 nm)
nv = NVCenter(NVspectrum.wavelength)

# Define objective and plate
obj = MPlanApo100x()
diamond_plate = Diamond(100.0)

# Compute polychromatic PSF at on-axis, z=0
psf_nv = PSF(0.0, 0.0, 0.0, obj, nv, diamond_plate)
println("NV-center weighted PSF intensity: ", psf_nv)
```

---

## Concluding Remarks

`VectorPSFs.jl` merges vectorial diffraction integrals, plate-dispersion modeling, and NV-center polychromatic weighting. The user can explore multiple incidence angles (`α`), different defocus values, and optimize over plate thickness or defocus to understand aberrations. The code is heavily based on numerical integration via `HCubature` and often uses `Optim` (Brent method) for 1D optimization tasks (finding best focus, searching for thickness).

**Key Exports**:

- **Objectives**: `MPlanApo100x`, `LMPLFLN100XBD`, `MPlanApo50x`, `M10x`, `UPLXAPO100X`
- **Plates**: `Diamond`, `FusedSilica`, `BorosilicateCrown`, `Sapphire`, `MagnesiumFluoride`, `CustomPlate`
- **PSF** function family
- **Strehl ratio** routines: `strehl(...)`
- **Thickness finder**: `maxtol_thick(...)`
- **NV**-center spline: `nv_spectrum_spline`, `NVCenter(...)`

Feel free to adjust tolerances (`rtol`, `atol`) as needed for speed vs. accuracy. If performing repeated calls (e.g., scanning `x, y, z` grids), you may wish to cache results or consider a higher-level approach to manage performance.

That covers the main functionality and usage of **VectorPSFs.jl**!