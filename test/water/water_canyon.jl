using Test
using UrbanTethysChloris.Water: water_canyon

FT = Float64

@testset "MATLAB" begin
    # Common parameters
    dth = 1.0
    row = 1000.0
    Zs = FT.([0, 10, 20, 50, 100, 150, 200, 300, 400, 600, 800, 1000, 1500, 2000])
    ms = length(Zs)-1

    # Create input NamedTuples
    MeteoData = (Rain=0.0,)
    Int_ittm = (
        IntGroundImp=0.0,
        IntGroundBare=0.0,
        IntGroundVegPlant=0.0,
        IntGroundVegGround=0.0,
        IntTree=0.0,
    )
    Owater_ittm = (
        OwGroundSoilImp=[NaN; NaN; fill(0.279610164940800, ms-2)],
        OwGroundSoilBare=fill(0.279610164940800, ms),
        OwGroundSoilVeg=fill(0.279610164940800, ms),
    )
    Runon_ittm = (RunonGroundTot=0.0,)
    Qinlat_ittm = (Qin_imp=zeros(FT, ms), Qin_bare=zeros(FT, ms), Qin_veg=zeros(FT, ms))
    ParSoilGround = (
        Pcla=0.2,
        Psan=0.4,
        Porg=0.025,
        Kfc=0.2,
        Phy=10000.0,
        SPAR=2,
        Kbot=NaN,
        Zs=Zs,
        In_max_imp=0.25,
        In_max_bare=10.0,
        In_max_underveg=10.0,
        Sp_In=0.2,
        Kimp=0.0,
    )
    ParInterceptionTree = (Sp_In=0.2,)
    ParCalculation = (dth=dth, row=row)
    ParVegGround = (
        LAI=2.5,
        SAI=1.0e-3,
        CASE_ROOT=1,
        ZR95=[250.0],
        ZR50=[NaN],
        ZRmax=[NaN],
        Rrootl=[3800.0],
        PsiL50=[-4.0],
        PsiX50=[-4.5],
    )
    ParVegTree = (
        LAI=5.0,
        SAI=0.2,
        CASE_ROOT=1,
        ZR95=[1000.0],
        ZR50=[NaN],
        ZRmax=[NaN],
        Rrootl=[4000.0],
        PsiL50=[-3.0],
        PsiX50=[-4.5],
        SPARTREE=2,
    )
    FractionsGround = (Per_runoff=0.9, fimp=0.455, fbare=0.0, fveg=0.545)
    geometry = (radius_tree=0.05,)
    ParTree = (trees=1,)
    Gemeotry_m = (Width_canyon=5.78,)
    Anthropogenic = (Waterf_canyonBare=0.0, Waterf_canyonVeg=0.0)

    # Evaporation/transpiration rates
    Egbare_Pond = 0.0
    Egbare_soil = 0.0
    Egimp_Pond = 0.0
    Egveg_In = 0.0
    Egveg_Pond = 0.0
    Egveg_soil = 2.431591473452235e-08
    Etree_In = 0.0
    TEgveg = 0.0
    TEtree = 0.0

    results = water_canyon(
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        Qinlat_ittm,
        Etree_In,
        Egveg_In,
        Egimp_Pond,
        Egbare_Pond,
        Egveg_Pond,
        Egbare_soil,
        Egveg_soil,
        TEgveg,
        TEtree,
        ParSoilGround,
        ParInterceptionTree,
        ParCalculation,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        geometry,
        ParTree,
        Gemeotry_m,
        Anthropogenic,
    )

    @test results.V_gimp ≈ [
        0,
        0,
        5.49247422004618,
        9.17015672378884,
        9.17363395996049,
        9.17430157201854,
        18.3486885448344,
        18.3486900465669,
        36.6973801302802,
        36.6973801304338,
        36.6973801304342,
        91.7434503260856,
        91.7434503260856,
    ] atol = 0.05
    @test results.Lk_gimp ≈ 0.017134579640331
    @test results.dV_dt_gimp ≈ -0.017077180764318 atol = 0.05
    @test results.Lk_gbare ≈ 0.017134579640331
    @test results.dV_dt_gbare ≈ -0.017134579640469
    @test results.Lk_gveg ≈ 0.017134579640331
    @test results.Runoff == 0
    @test results.Etot ≈ 4.770782470913287e-05
    @test results.DeepGLk ≈ 0.017134579640331
    @test results.StorageTot ≈ -0.017182287465060
end
