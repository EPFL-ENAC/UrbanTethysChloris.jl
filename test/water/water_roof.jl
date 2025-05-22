using Test
using UrbanTethysChloris.Water: water_roof

FT = Float64

@testset "MATLAB - Zurich" begin
    Eroof_imp = -2.182498032323756e-10
    Eroof_veg = 0.0
    Eroof_ground = 0.0
    Eroof_soil = 5.746195678496601e-10
    TEroof_veg = 0.0

    MeteoData = (Rain=0.0,)
    Int_ittm = (IntRoofImp=0.0, IntRoofVegPlant=0.0, IntRoofVegGround=0.0)
    Owater_ittm = (OwRoofSoilVeg=fill(0.279610164940800, 4),)
    Runon_ittm = (RunonRoofTot=0.0,)
    FractionsRoof = (Per_runoff=1.0, fveg=0.5, fimp=0.5)
    ParSoilRoof = (
        In_max_imp=0.25,
        In_max_ground=10.0,
        Kimp=0.0,
        Sp_In=0.2,
        Zs=[0.0, 10.0, 20.0, 50.0, 100.0],
        Pcla=0.2,
        Psan=0.4,
        Porg=0.025,
        Kfc=0.2,
        Phy=10000.0,
        SPAR=2,
        Kbot=NaN,
    )
    ParCalculation = (dth=1.0, row=1000.0)

    ParVegRoof = (
        LAI=2.5,
        SAI=1.0e-3,
        Rrootl=[3800.0],
        PsiL50=[-4.0],
        PsiX50=[-4.5],
        CASE_ROOT=1,
        ZR95=[70.0],
        ZR50=[NaN],
        ZRmax=[NaN],
    )
    Anthropogenic = (Waterf_roof=0.0,)

    q_runon_imp, In_imp, dIn_imp_dt, Lk_imp, q_runon_veg, In_veg, dIn_veg_dt, q_runon_ground, In_ground, dIn_ground_dt, dV_dt, f_ground, V, O, OS, Lk, Psi_s, Exwat, Rd, TEroof_veg, Eroof_soil, Runoff, Runon, WBalance_In_imp, WBalance_In_veg, WBalance_In_ground, WBalance_soil, WBalance_imp_tot, WBalance_veg_tot, WBalance_tot = water_roof(
        Eroof_imp,
        Eroof_veg,
        Eroof_ground,
        Eroof_soil,
        TEroof_veg,
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        FractionsRoof,
        ParSoilRoof,
        ParCalculation,
        ParVegRoof,
        Anthropogenic,
    )

    @test q_runon_imp ≈ 0.0
    @test In_imp ≈ 7.856992916365522e-07
    @test dIn_imp_dt ≈ 7.856992916365522e-07
    @test Lk_imp ≈ 0.0
    @test q_runon_veg ≈ 0.0
    @test In_veg ≈ 0.0
    @test dIn_veg_dt ≈ 0.0
    @test q_runon_ground ≈ 0.0
    @test In_ground ≈ 0.0
    @test dIn_ground_dt ≈ 0.0
    @test dV_dt ≈ -0.017136648270778
    @test f_ground ≈ 0.0
    @test V ≈ [1.829564658402304, 1.830898843370262, 5.498728268305371, 9.172361646868405] atol =
        0.05
    @test O ≈ [0.279079730128859, 0.279213148625655, 0.279414206565475, 0.279570497225997] atol =
        0.05
    @test OS ≈ 0.279079730128859 atol = 0.05
    @test Lk ≈ 0.017134579640331 atol = 0.05
    @test Psi_s ≈ [-0.033225482352765] atol = 0.05
    # @test Exwat ≈ [48023.64755288577, 31490.63736321513, 43007.01370915098, 9531.74789365364] atol = odetol
    @test Rd == 0
    @test TEroof_veg == 0
    @test Eroof_soil ≈ 5.746195678496601e-10
    @test Runoff ≈ 0.0
    @test Runon ≈ 0.0
    @test WBalance_In_imp == 0
    @test WBalance_In_veg == 0
    @test WBalance_In_ground == 0
    @test WBalance_soil ≈ 0.0 atol = 1e-12
    @test WBalance_imp_tot == 0
    @test WBalance_veg_tot ≈ 0.0 atol = 1e-12
    @test WBalance_tot ≈ 0.0 atol = 1e-12
end
