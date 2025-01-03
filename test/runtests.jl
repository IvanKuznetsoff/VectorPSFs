using VectorPSFs
using Test
using Statistics
import PSFModels as M

@testset "VectorPSFs.jl" begin
    # Write your tests here.
    FWHM = 1.434638
    model = M.airydisk(x = 0, y=0, fwhm = FWHM)
    p0 = PSFParams(0.7, M10x())
    p1 = PSFParams(0.7, M10x(), Diamond(0.0))
    p2 = PSFParams(0.7, M10x(), Diamond(0.0), 0.0)

    xs = range(-3,3,1001)
    ys = [PSF(x,0,0,p0; atol = 1e-6, rtol = 1e-5) for x in xs] ./ PSF(0,0,0,p0; atol = 1e-6, rtol = 1e-5) 
    @test sqrt(mean(@. abs2(ys .- model.(xs, 0)))) < 5e-4

    ys .= [PSF(x,0,0,p1; atol = 1e-6, rtol = 1e-5) for x in xs] ./ PSF(0,0,0,p1; atol = 1e-6, rtol = 1e-5) 
    @test sqrt(mean(@. abs2(ys .- model.(xs, 0)))) < 5e-4
    
    ys .= [PSF(x,0,0,p2; atol = 1e-6, rtol = 1e-5) for x in xs] ./ PSF(0,0,0,p2; atol = 1e-6, rtol = 1e-5) 
    @test sqrt(mean(@. abs2(ys .- model.(xs, 0)))) < 5e-4

    @test round(maxtol_thick(0.7, MPlanApo100x(), Diamond)) ≈ 41
    @test round(maxtol_thick(0.7, MPlanApo100x(), Diamond, 0.0)) ≈ 41
    # @test round(maxtol_thick(0.7, MPlanApo100x(), Diamond, asin(100/2000))) ≈ 31

    nv = NVCenter([0.7], [1.])
    obj = M10x()
    
    ys = [PSF(x, 0, 0, obj, nv, Diamond(0.0), 0.0; atol = 1e-6, rtol = 1e-5) for x in xs] ./ PSF(0, 0, 0, obj, nv, Diamond(0.0), 0.0; atol = 1e-6, rtol = 1e-5) 
    @test sqrt(mean(@. abs2(ys .- model.(xs, 0)))) < 5e-4

end
