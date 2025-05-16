using Test
using UrbanTethysChloris.Soil: infiltration

@testset "MATLAB" begin
    alpVG = -9.634329919190001e-05
    cosalp = 1.0
    Dz = 10.0
    Ks_Zs = 15.475656399419266
    L = 0.186873685240403
    nVG = 1.186873685240403
    O = 0.279610164940800
    O33 = 0.279610164940800
    Ohy = 0.096123264288629
    Osat = 0.459478244940800
    Pe = 3.912121076244333
    Pond = 0.0
    SPAR = 2
    WIS = 0.0

    f, fpot = infiltration(
        Osat, Ohy, L, alpVG, nVG, Pe, Ks_Zs, O33, SPAR, O, Dz, WIS, cosalp, Pond
    )

    # Values should be positive and fpot should be >= f
    @test f == 0
    @test fpot ≈ 2304.649699209915
end
