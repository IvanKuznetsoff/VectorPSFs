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
using VectorPsfs

obj = MPlanApo50x()         # e.g. 50x objective
λ   = 0.632                 # 632 nm = 0.632 µm
```

To obtain a PSF map on a plane:
```julia
xs = range(-1,1,201)
ys = range(-1,1,201)
psf_map_xy = [PSF(x, y, 0.0, λ, obj) for x in xs, y in ys]
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
using VectorPsfs, Optim

obj = MPlanApo100x()
λ   = 0.7         # 700 nm
s = Strehl(50.0, λ, obj, Diamond; zrange=[0.0, 5.0])  # Diamond thickness=50 µm
println("Strehl ratio = ", s)
```

## NV-Center Weighted PSF

We also offer a polychromatic approach for NV emission:  
```julia
PSF(x, y, z, obj::Objective, nv::NVCenter; ...)
```
where `NVCenter` is a struct storing smoothed spectral weights. This integrates over multiple wavelengths.  

**Example**:
```julia
using VectorPsfs

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
```
