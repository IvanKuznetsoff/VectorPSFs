# VectorPSFs

[![Build Status](https://github.com/IvanKuznetsoff/VectorPSFs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/IvanKuznetsoff/VectorPSFs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://app.travis-ci.com/IvanKuznetsoff/VectorPSFs.jl.svg?branch=main)](https://app.travis-ci.com/IvanKuznetsoff/VectorPSFs.jl)
[![Coverage](https://codecov.io/gh/IvanKuznetsoff/VectorPSFs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/IvanKuznetsoff/VectorPSFs.jl)
[![Coverage](https://coveralls.io/repos/github/IvanKuznetsoff/VectorPSFs.jl/badge.svg?branch=main)](https://coveralls.io/github/IvanKuznetsoff/VectorPSFs.jl?branch=main)

## Overview

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

For installation and setup instructions, see [Document](https://ivankuznetsoff.github.io/VectorPSFs.jl/).
