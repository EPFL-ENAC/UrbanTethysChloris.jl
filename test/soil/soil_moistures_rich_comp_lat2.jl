using Test
using UrbanTethysChloris.Soil: soil_moistures_rich_comp_lat2

FT = Float64

@testset "MATLAB" begin
    alpVG = -9.634329919190001e-05
    C1 = 1.0
    C2 = 1.0
    dz = 10.0
    f1 = 0.1
    f2 = 0.23
    Ks = 15.475656399419266
    L = 0.186873685240403
    nVG = 1.186873685240403
    O33 = 0.279610164940800
    Ohy = 0.096123264288629
    Osat = 0.459478244940800
    Pe = 3.912121076244333
    SPAR = 2
    Vlat = [1.806497419003835, 1.806335223222684]
    Wcan = 5.78

    #
    @testset "SPAR 2" begin
        dVlat = soil_moistures_rich_comp_lat2(
            Vlat, dz, SPAR, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, C1, C2, f1, f2, Wcan
        )

        @test dVlat ≈ [-4.292066415728490e-06, 1.866115832925431e-06] atol = 1e-7
    end

    @testset "SPAR 1" begin
        dVlat = soil_moistures_rich_comp_lat2(
            Vlat, dz, 1, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, C1, C2, f1, f2, Wcan
        )

        @test dVlat ≈ [-2.073296129980159e-06, 9.014330999913737e-07] atol = 1e-8
    end
end
