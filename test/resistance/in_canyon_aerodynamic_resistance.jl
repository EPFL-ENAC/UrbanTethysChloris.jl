using Test
using UrbanTethysChloris.Resistance: in_canyon_aerodynamic_resistance

FT = Float64

@testset "MATLAB" begin
    dcan = 4.902710934818196
    ea = 763.6587960785345
    hcan = 6.5
    hcan_max = NaN
    Pre = 94872.32105699336
    RoughnessParameter = :MacD
    Ta = 3.9
    Ts = 3.9
    uatm = 0.4
    Zatm = 25.0
    zom_und = 0.01845
    zomcan = 0.364177306663220
    Zp1 = 4.711
    Zp2 = 2.0
    Zp3 = 2.0
    Zref_und = 1.5

    rap_can, rap_Zp1, rap_Zp1_In, rap_Zp2, rap_Zp2_In, rap_Zp3, rap_Zp3_In, u_Hcan, u_Zp1, u_Zp2, u_Zp3, uref_und, alpha = in_canyon_aerodynamic_resistance(
        uatm,
        Zatm,
        Ta,
        Ts,
        hcan_max,
        hcan,
        dcan,
        zomcan,
        Zref_und,
        zom_und,
        Zp1,
        Zp2,
        Zp3,
        Pre,
        ea,
        RoughnessParameter,
    )

    @test rap_can ≈ 1188.766421761054
    @test rap_Zp1 ≈ 1166.24405858974
    @test rap_Zp1_In ≈ 22.522363171314055
    @test rap_Zp2 ≈ 1046.184893171181
    @test rap_Zp2_In ≈ 142.5815285898730
    @test rap_Zp3 ≈ 1046.184893171181
    @test rap_Zp3_In ≈ 142.5815285898730
    @test u_Hcan ≈ 0.147447828918469
    @test u_Zp1 ≈ 0.133882970988947
    @test u_Zp2 ≈ 0.115667541369742
    @test u_Zp3 ≈ 0.115667541369742
    @test uref_und ≈ 0.112589374601137
    @test alpha ≈ 0.350645182969221
end
