using Test
using UrbanTethysChloris.Soil: soil_parameters2

FT = Float64

@testset "MATLAB" begin
    Osat = [0.459478244940800]
    L = [0.186873685240403]
    Pe = [3.912121076244333]
    Ks = [15.475656399419266]
    O33 = [0.279610164940800]
    nVG = [1.5]
    alpVG = [0.02]
    Kfc = 1.0
    Pss = 100.0
    Pwp = 15000.0
    Phy = 10000.0

    @testset "SPAR 2" begin
        Ofc, Oss, Owp, Ohy = soil_parameters2(
            Osat, L, Pe, Ks, O33, nVG, alpVG, Kfc, Pss, Pwp, Phy
        )

        @test all(Ofc .≈ 0.376222329280883)
        @test all(Oss .≈ 0.227287672446908)
        @test all(Owp .≈ 0.089109024271719)
        @test all(Ohy .≈ 0.096123264288629)
    end

    @testset "SPAR 1" begin
        Ohy = [0.042]
        Ofc, Oss, Owp, Ohy = soil_parameters2(
            Osat, L, Pe, Ks, O33, nVG, alpVG, Kfc, Pss, Pwp, Phy, Ohy, 1
        )

        @test all(Ofc .== 0)
        @test all(Oss .≈ 0.071235035891368)
        @test all(Owp .≈ 0.044387303839774)
        @test all(Ohy .== 0)
    end
end
