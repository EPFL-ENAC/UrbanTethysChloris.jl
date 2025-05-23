using Test
using UrbanTethysChloris.Resistance: wind_profile_roof

FT = Float64

@testset "MATLAB" begin
    disp_h = 0.1
    Hcan = 6.5
    hveg = 0.15
    uatm = 0.4
    Zatm = 25.0
    zom = 0.01845
    Zp = 2.0

    u_Zp, u_Hveg = wind_profile_roof(Hcan, Zatm, uatm, zom, disp_h, hveg, Zp)

    @test u_Zp ≈ 4.122087360758734
    @test u_Hveg ≈ 2.302585092994046
end
