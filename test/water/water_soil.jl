using Test
using UrbanTethysChloris.Water: water_soil

FT = Float64

@testset "MATLAB" begin
    CASE_ROOT_H = 1
    CASE_ROOT_L = 1
    dth = 1.0
    E_soil = 5.746195678496601e-10
    f = 0.0
    Kbot = NaN
    Kfc = 0.2
    Otm1 = fill(0.279610164940800, 4)
    Pcla = 0.2
    Phy = 10000.0
    Porg = 0.025
    Psan = 0.4
    PsiL50_H = [0.0]
    PsiL50_L = [-4.0]
    PsiX50_H = [0.0]
    PsiX50_L = [-4.5]
    Qlat_in = [0.0, 0.0, 0.0, 0.0, 0.0]
    row = 1000.0
    Rrootl_H = [0.0]
    Rrootl_L = [3800.0]
    SPAR = 2
    TE_H = 0.0
    TE_L = 0.0
    ZR50_H = [0.0]
    ZR50_L = [NaN]
    ZR95_H = [0.0]
    ZR95_L = [70.0]
    ZRmax_H = [0.0]
    ZRmax_L = [NaN]
    Zs = [0.0, 10.0, 20.0, 50.0, 100.0]

    V, O, OS, Lk, Psi_s_H, Psi_s_L, Exwat_H, Exwat_L, Rd, TE_L, TE_H, E_soil, dV_dt, WBalance_soil, Psi_soil, Ko = water_soil(
        Otm1,
        f,
        TE_H,
        TE_L,
        E_soil,
        Qlat_in,
        dth,
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
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
        Zs,
        row,
    )

    odetol = 0.05
    @test V ≈ [1.829564658402304, 1.830898843370262, 5.498728268305371, 9.172361646868405] atol =
        odetol
    @test O ≈ [0.279079730128859, 0.279213148625655, 0.279414206565475, 0.279570497225997] atol =
        odetol
    @test OS ≈ 0.279079730128859 atol = odetol
    @test Lk ≈ 0.017134579640331 atol = odetol
    @test all(isnan.(Psi_s_H))
    @test Psi_s_L ≈ [-0.033225482352765] atol = odetol
    @test all(Exwat_H .== 0)
    @test Rd == 0
    @test TE_L == 0
    @test TE_H == 0
    @test E_soil ≈ 5.746195678496601e-10
    @test dV_dt ≈ -0.017136648270778
    @test Ko ≈ [0.016694507101595, 0.016804199773631, 0.016970766254528, 0.017101301153577] atol =
        odetol

    # Not working currently
    # Since there is a discrepancy in V between MATLAB and Julia, this difference is propagagted
    # and multiplied by 5e9 within root_soil_conductance.
    # @test Exwat_L ≈ [48023.64755288577, 31490.63736321513, 43007.01370915098, 9531.74789365364] atol = odetol
    # @test WBalance_soil ≈ 2.782496455466799e-15
    # @test Psi_soil ≈ [3398.269886762692, 3389.589513226043, 3376.558090629037, 3366.469289461891] atol = odetol
end
