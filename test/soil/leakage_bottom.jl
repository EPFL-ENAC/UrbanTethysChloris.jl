using Test
using UrbanTethysChloris.Soil: leakage_bottom

FT = Float64

@testset "MATLAB" begin
    ms = 1
    Ks_Zs = [15.475656399419266]
    L = [0.186873685240403]
    nVG = [1.186873685240403]
    O = [0.279610164940800]
    Ohy = [0.096123264288629]
    Osat = [0.459478244940800]

    # Test with no bottom conductivity
    @testset "No bottom conductivity" begin
        Kbot = NaN
        Lk = leakage_bottom(O, Ks_Zs, Osat, Ohy, L, nVG, Kbot, ms, 1)
        @test Lk ≈ 4.690660454201376e-05

        Lk = leakage_bottom(O, Ks_Zs, Osat, Ohy, L, nVG, Kbot, ms, 2)
        @test Lk ≈ 0.017134579640331
    end

    @testset "Test with bottom conductivity" begin
        Kbot = 0.1
        @test leakage_bottom(O, Ks_Zs, Osat, Ohy, L, nVG, Kbot, ms, 2) == 0.0

        O[ms] = Osat[ms]
        @test leakage_bottom(O, Ks_Zs, Osat, Ohy, L, nVG, Kbot, ms, 2) == Kbot
    end
end
