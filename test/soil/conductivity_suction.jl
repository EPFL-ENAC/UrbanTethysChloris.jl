using Test
using UrbanTethysChloris.Soil: conductivity_suction

@testset "MATLAB" begin
    alpVG = -9.634329919190001e-05
    Ks = 15.475656399419266
    L = 0.186873685240403
    nVG = 1.186873685240403
    O = 0.279610164940800
    O33 = 0.279610164940800
    Ohy = 0.096123264288629
    Osat = 0.459478244940800
    Pe = 3.912121076244333

    @testset "SPAR 1" begin
        Ko, Po = conductivity_suction(1, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, O)

        @test Ko ≈ 4.690660454201376e-05
        @test Po ≈ 397389.8386717513
    end

    @testset "SPAR 2" begin
        Ko, Po = conductivity_suction(2, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, O)

        @test Ko ≈ 0.017134579640331
        @test Po ≈ 3363.914373088685
    end
end
