using Test
using UrbanTethysChloris.Resistance: wind_profile_point_output

FT = Float64

@testset "MATLAB" begin
    Zp = 1.1

    Gemeotry_m = (
        Height_canyon=6.5,
        Width_canyon=5.78,
        Width_roof=4.73,
        Height_tree=4.711,
        Radius_tree=0.289,
        Hcan_max=NaN,
        Hcan_std=NaN,
    )

    ParVegTree = (Kopt=0.5, LAI=5.0)

    ParTree = (trees=true,)

    MeteoData = (Zatm=25.0, Uatm=0.4)

    FractionsGround = (fveg=0.455, fbare=0.0, fimp=0.55)

    ParVegGround = (hc=0.15,)

    u_Zp = wind_profile_point_output(
        Zp, Gemeotry_m, ParVegTree, ParTree, MeteoData, FractionsGround, ParVegGround
    )

    @test u_Zp ≈ 0.104649649571975
end
