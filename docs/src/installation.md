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

Once installed, you can `using VectorPsfs` in your Julia session:
```julia
using VectorPsfs
```

Check your environment’s compatibility with **Julia ≥ 1.10**.
```

