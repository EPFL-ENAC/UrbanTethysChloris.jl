using Test
using UrbanTethysChloris.Resistance: wind_profile_canyon

FT = Float64

@testset "MATLAB" begin
    Hcan = 6.5
    Hcan_max = NaN
    Hcan_std = NaN
    Htree = 4.711
    Kopt = 0.5
    LAI_t = 5.0
    R_tree = 0.289
    trees = true
    uatm = 0.4
    Wcan = 5.78
    Wroof = 4.73
    Zatm = 25.0
    zom_und = 0.01
    Zp = 6.5
    Zref_und = 1.5

    dcan, zomcan, u_Hcan_max, u_Zp, w_Zp, alpha, RoughnessParameter = wind_profile_canyon(
        Hcan,
        Htree,
        R_tree,
        Wcan,
        Wroof,
        Kopt,
        LAI_t,
        Zatm,
        uatm,
        Zp,
        trees,
        Zref_und,
        zom_und,
        Hcan_max,
        Hcan_std,
    )

    @test dcan ≈ 4.902710934818196
    @test zomcan ≈ 0.364177306663220
    @test u_Hcan_max ≈ 0.147447828918469
    @test u_Zp ≈ 0.147447828918469
    @test w_Zp == 0
    @test alpha ≈ 0.350645182969221
    @test RoughnessParameter == :MacD
end
