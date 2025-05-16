using Test
using UrbanTethysChloris.Soil: soil_moistures_rich_comp_lat3

FT = Float64

@testset "MATLAB" begin
    alpVG = -9.634329919190001e-05
    Cimp = 1.0
    Cbare = 1.0
    Cveg = 1.0
    dz = 30.0
    fbare = 0.1
    fimp = 0.67
    fveg = 0.23
    Ks = 15.475656399419266
    L = 0.186873685240403
    nVG = 1.186873685240403
    O33 = 0.279610164940800
    Ohy = 0.096123264288629
    Osat = 0.459478244940800
    Pe = 3.912121076244333
    Vlat = [5.492435228207039, 5.469208727510397, 5.469009151102071]
    Wcan = 5.8

    @testset "SPAR 2" begin
        dVlat = soil_moistures_rich_comp_lat3(
            Vlat,
            dz,
            2,
            Ks,
            Osat,
            Ohy,
            L,
            Pe,
            O33,
            alpVG,
            nVG,
            Cimp,
            Cbare,
            Cveg,
            fimp,
            fbare,
            fveg,
            Wcan,
        )

        @test dVlat ≈ [-0.0001945704451840972, 0.0006435312679622614, 0.0002869959629439957]
    end

    @testset "SPAR 1" begin
        dVlat = soil_moistures_rich_comp_lat3(
            Vlat,
            dz,
            1,
            Ks,
            Osat,
            Ohy,
            L,
            Pe,
            O33,
            alpVG,
            nVG,
            Cimp,
            Cbare,
            Cveg,
            fimp,
            fbare,
            fveg,
            Wcan,
        )

        @test dVlat ≈ [-9.628807920288298e-05, 0.0003184861276097109, 0.0001420191317606978]
    end
end
