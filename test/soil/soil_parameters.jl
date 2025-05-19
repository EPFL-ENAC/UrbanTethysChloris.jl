using Test
using UrbanTethysChloris.Soil: soil_parameters

FT = Float64

@testset "MATLAB" begin
    Psan = 0.4
    Pcla = 0.2
    Porg = 0.025

    Osat, L, Pe, Ks, O33, rsd, lan_dry, lan_s, cv_s, K_usle = soil_parameters(
        Psan, Pcla, Porg
    )

    @test Osat ≈ 0.459478244940800
    @test L ≈ 0.186873685240403
    @test Pe ≈ 3.912121076244333
    @test Ks ≈ 15.475656399419266
    @test O33 ≈ 0.279610164940800
    @test rsd ≈ 1432.382650906880
    @test lan_dry ≈ 0.191773412109454
    @test lan_s ≈ 6.537988142956436
    @test cv_s ≈ 2.226788982559455e6
    @test K_usle ≈ 4.339377433037157e-5

    @test_throws ArgumentError soil_parameters(1.0, 1.0, 1.0)
    @test_throws ArgumentError soil_parameters(-1.0, 1.0, 1.0)
end
