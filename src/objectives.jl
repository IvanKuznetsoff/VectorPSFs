# objectives.jl

# Struct definitions and constructors for microscope objectives

struct Objective
    f::Float64    # Focal length
    NA::Float64   # Numerical aperture
    n::Float64    # Refractive index in vacuum (often 1.0 for air)
end


# Mitsutoyo MPlanApo100x returns an Objective with 
# f=2000, NA=0.7, n=1.0

function MPlanApo100x()
    return Objective(2000, 0.7, 1.0)
end


# Olympus LMPLFLN100XBD returns an Objective with
# f=1800, NA=0.8, n=1.0

function LMPLFLN100XBD()
    return Objective(1800, 0.8, 1.0)
end


# Mitsutoyo MPlanApo50x returns an Objective with
# f=4000, NA=0.55, n=1.0

function MPlanApo50x()
    return Objective(4000, 0.55, 1.0)
end


# Newport M10x returns an Objective with 
# f=16500, NA=0.25, n=1.0

function M10x()
    return Objective(16500, 0.25, 1.0)
end

# Olympus UPLXAPO100X returns an Objective with
# f=1800, NA=1.45, n=1.45

function UPLXAPO100X()
    return Objective(1800, 1.45, 1.515)
end