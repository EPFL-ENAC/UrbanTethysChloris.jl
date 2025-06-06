using Test
using UrbanTethysChloris.Resistance: soil_resistance

FT = Float64

@testset "MATLAB" begin
    alpVG = -9.634329919190001e-05
    ea = 763.6587960785345
    Ks = 15.475656399419266
    L = 0.186873685240403
    nVG = 1.186873685240403
    O = 0.279610164940800
    O33 = 0.279610164940800
    Ohy = 0.096123264288629
    Osat = 0.459478244940800
    Pe = 3.912121076244333
    Pre = 948.7232105699336
    q_runon = 0.0
    SPAR = 2
    Ts = 3.899999999999977
    Ws = 0.268262529232655

    r_soil, b_soil, alp_soil = soil_resistance(
        Ts, Pre, Ws, ea, q_runon, O, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, SPAR
    )
    @test r_soil ≈ 233.318446945321
    @test b_soil ≈ 1
    @test alp_soil ≈ 0.999741935636661
end
