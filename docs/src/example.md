# Examples

`PSF(x, y, z, ...)` is defined in multiple overloads:

- **No plate, direct calculation** (simple, air):
  ```julia
  PSF(x, y, z, λ, NA; rtol=1e-3, atol=1e-4)
  ```
- **No plate, with objective specification** (aberration-free, immersion):
  ```julia
  PSF(x, y, z, λ, obj::Objective; rtol=1e-3, atol=1e-4)
  ```
- **Normal-incidence plate, with objective**:
  ```julia
  PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate; ...)
  ```
- **Tilted-incidence plate, with objective**:
  ```julia
  PSF(x, y, z, λ, obj::Objective, plate::PlaneParallelPlate, α; ...)
  ```

`(x, y, z)` are spatial coordinates (in μm), `λ` is wavelength in micrometers, `NA` is the numerical aperture, `α` is the tilt angle in radians.

The kwargs `rtol`/`atol` control the tolerance in the numerical integration (via `HCubature`).

**Example** (Normal incidence, plate, objective):
```julia
using VectorPSFs

obj = MPlanApo50x()           # e.g. 50x objective
plate = Diamond(50.0)         # 50 μm diamond
λ = 0.632                     # 632 nm = 0.632 μm
psf_val = PSF(0.0, 0.0, 0.0, λ, obj, plate)  
println("PSF on-axis intensity: ", psf_val)
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