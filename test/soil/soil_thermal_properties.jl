using Test
using UrbanTethysChloris.Soil: soil_thermal_properties

FT = Float64

@testset "MATLAB" begin
    O = fill(0.279610164940800, 4)
    Ohy = fill(0.096123264288629, 4)
    Osat = fill(0.459478244940800, 4)
    Tdp = fill(3.899999999999977, 4)
    cv_s = fill(2.213666666666667e+06, 4)
    lan_dry = fill(0.198582784280192, 4)
    lan_s = fill(6.839999999999999, 4)
    rsd = fill(1432.38265090688, 4)

    lanS, cv_Soil, CTt = soil_thermal_properties(
        Tdp, rsd, lan_dry, lan_s, cv_s, Osat, Ohy, O
    )

    @test all(lanS .≈ 1.206526914774966)
    @test all(cv_Soil .≈ 2.366983142224905e+06)
    @test all(CTt .≈ 7.136445789505452e-06)
end
