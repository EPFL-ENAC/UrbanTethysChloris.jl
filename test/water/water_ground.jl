using Test
using UrbanTethysChloris.Water: water_ground

FT = Float64

@testset "MATLAB" begin
    CASE_ROOT_H = 1
    CASE_ROOT_L = 1
    dth = 1.0
    E_ground = 0.0
    In_ground_tm1 = 0.0
    In_max_ground = 10.0
    Kbot = 1.0
    Kfc = 0.2
    Otm1 = zeros(FT, 4)
    Pcla = 0.2
    Phy = 10000.0
    Porg = 0.025
    Psan = 0.4
    q_runon_veg = 0.0
    row = 1000.0
    Runon_tm1 = 0.0
    SPAR = 2
    ZR50_H = [0.0]
    ZR50_L = [NaN]
    ZR95_H = [0.0]
    ZR95_L = [95.0]
    ZRmax_H = [0.0]
    ZRmax_L = [NaN]
    Zs = [0.0, 10.0, 20.0, 50.0, 100.0]

    q_runon_ground, In_ground, dIn_ground_dt, f_ground, WBalance_In_ground = water_ground(
        q_runon_veg,
        Runon_tm1,
        E_ground,
        Otm1,
        In_ground_tm1,
        In_max_ground,
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT_H,
        CASE_ROOT_L,
        ZR95_H,
        ZR95_L,
        ZR50_H,
        ZR50_L,
        ZRmax_H,
        ZRmax_L,
        Zs,
        dth,
        row,
    )

    @test q_runon_ground ≈ 0.0
    @test In_ground ≈ 0.0
    @test dIn_ground_dt ≈ 0.0
    @test f_ground ≈ 0.0
    @test WBalance_In_ground ≈ 0.0
end
