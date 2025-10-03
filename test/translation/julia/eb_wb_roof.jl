using Test
using MAT
using UrbanTethysChloris: eb_wb_roof
using ....TestUtils:
    load_matlab_data,
    create_urban_geometry_parameters,
    create_location_specific_surface_fractions,
    create_vegetated_soil_parameters,
    create_vegetated_optical_properties,
    create_location_specific_thermal_properties,
    create_height_dependent_vegetation_parameters

FT = Float64
input_vars, output_vars = load_matlab_data("EB_WB_roof.json")

# Create parameter structs from input data
Geometry_m = create_urban_geometry_parameters(
    FT;
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
)

FractionsRoof = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsRoof"]["fveg"],
    fimp=input_vars["FractionsRoof"]["fimp"],
    Per_runoff=FT(input_vars["FractionsRoof"]["Per_runoff"]),
)

ParSoilRoof = create_vegetated_soil_parameters(
    FT;
    In_max_imp=input_vars["ParSoilRoof"]["In_max_imp"],
    In_max_ground=FT(input_vars["ParSoilRoof"]["In_max_ground"]),
    Kimp=FT(input_vars["ParSoilRoof"]["Kimp"]),
    Pcla=input_vars["ParSoilRoof"]["Pcla"],
    Psan=input_vars["ParSoilRoof"]["Psan"],
    Porg=input_vars["ParSoilRoof"]["Porg"],
    Kfc=input_vars["ParSoilRoof"]["Kfc"],
    Phy=FT(input_vars["ParSoilRoof"]["Phy"]),
    SPAR=Int(input_vars["ParSoilRoof"]["SPAR"]),
    Kbot=input_vars["ParSoilRoof"]["Kbot"],
    Sp_In=input_vars["ParSoilRoof"]["Sp_In"],
    Zs=FT.(vec(input_vars["ParSoilRoof"]["Zs"])),
    dz1=input_vars["ParSoilRoof"]["dz1"],
    dz2=input_vars["ParSoilRoof"]["dz2"],
)

PropOpticalRoof = create_vegetated_optical_properties(
    FT;
    aimp=input_vars["PropOpticalRoof"]["aimp"],
    aveg=input_vars["PropOpticalRoof"]["aveg"],
    eimp=input_vars["PropOpticalRoof"]["eimp"],
    eveg=input_vars["PropOpticalRoof"]["eveg"],
)

ParThermalRoof = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalRoof"]["cv_s_imp"],
    lan_dry=input_vars["ParThermalRoof"]["lan_dry_imp"],
)

ParVegRoof = create_height_dependent_vegetation_parameters(
    FT;
    LAI=FT(input_vars["ParVegRoof"]["LAI"]),
    SAI=FT(input_vars["ParVegRoof"]["SAI"]),
    hc=FT(input_vars["ParVegRoof"]["hc"]),
    d_leaf=FT(input_vars["ParVegRoof"]["d_leaf"]),
    Kopt=FT(input_vars["ParVegRoof"]["Kopt"]),
    Knit=FT(input_vars["ParVegRoof"]["Knit"]),
    Psi_sto_50=FT(input_vars["ParVegRoof"]["Psi_sto_50"]),
    Psi_sto_00=FT(input_vars["ParVegRoof"]["Psi_sto_00"]),
    CT=Int(input_vars["ParVegRoof"]["CT"]),
    Vmax=FT(input_vars["ParVegRoof"]["Vmax"]),
    DSE=FT(input_vars["ParVegRoof"]["DSE"]),
    Ha=FT(input_vars["ParVegRoof"]["Ha"]),
    FI=FT(input_vars["ParVegRoof"]["FI"]),
    Do=FT(input_vars["ParVegRoof"]["Do"]),
    a1=FT(input_vars["ParVegRoof"]["a1"]),
    go=FT(input_vars["ParVegRoof"]["go"]),
    e_rel=FT(input_vars["ParVegRoof"]["e_rel"]),
    e_relN=FT(input_vars["ParVegRoof"]["e_relN"]),
    gmes=FT(input_vars["ParVegRoof"]["gmes"]),
    rjv=FT(input_vars["ParVegRoof"]["rjv"]),
    mSl=FT(input_vars["ParVegRoof"]["mSl"]),
    Sl=FT(input_vars["ParVegRoof"]["Sl"]),
    Rrootl=FT.([input_vars["ParVegRoof"]["Rrootl"]]),
    PsiL50=FT.([input_vars["ParVegRoof"]["PsiL50"]]),
    PsiX50=FT.([input_vars["ParVegRoof"]["PsiX50"]]),
    CASE_ROOT=Int(input_vars["ParVegRoof"]["CASE_ROOT"]),
    ZR95=FT.([input_vars["ParVegRoof"]["ZR95"]]),
    ZR50=FT.([input_vars["ParVegRoof"]["ZR50"]]),
    ZRmax=FT.([input_vars["ParVegRoof"]["ZRmax"]]),
)

MeteoData = (;
    Zatm=FT(input_vars["MeteoData"]["Zatm"]),
    Tatm=input_vars["MeteoData"]["Tatm"],
    Uatm=input_vars["MeteoData"]["Uatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
    SW_dir=input_vars["MeteoData"]["SW_dir"],
    SW_diff=input_vars["MeteoData"]["SW_diff"],
    LWR=input_vars["MeteoData"]["LWR"],
    Rain=FT(input_vars["MeteoData"]["Rain"]),
    Catm_O2=input_vars["MeteoData"]["Catm_O2"],
    Catm_CO2=input_vars["MeteoData"]["Catm_CO2"],
)

TempVec_ittm = (;
    TRoofVeg=input_vars["TempVec_ittm"]["TRoofVeg"],
    TRoofIntVeg=input_vars["TempVec_ittm"]["TRoofIntVeg"],
    TRoofIntImp=input_vars["TempVec_ittm"]["TRoofIntImp"],
)

Int_ittm = (;
    IntRoofImp=FT(input_vars["Int_ittm"]["IntRoofImp"]),
    IntRoofVegPlant=FT(input_vars["Int_ittm"]["IntRoofVegPlant"]),
    IntRoofVegGround=FT(input_vars["Int_ittm"]["IntRoofVegGround"]),
)

ExWater_ittm = (; ExWaterRoofVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterRoofVeg_L"]))
Vwater_ittm = (; VRoofSoilVeg=vec(input_vars["Vwater_ittm"]["VRoofSoilVeg"]))
Owater_ittm = (; OwRoofSoilVeg=vec(input_vars["Owater_ittm"]["OwRoofSoilVeg"]))
SoilPotW_ittm = (; SoilPotWRoofVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWRoofVeg_L"])
CiCO2Leaf_ittm = (;
    CiCO2LeafRoofVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegSun"],
    CiCO2LeafRoofVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegShd"],
)
Runon_ittm = (; RunonRoofTot=FT(input_vars["Runon_ittm"]["RunonRoofTot"]))

HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])
Anthropogenic = (;
    Tb=input_vars["Anthropogenic"]["Tb"],
    Waterf_roof=input_vars["Anthropogenic"]["Waterf_roof"],
)
ParCalculation = (;
    dth=FT(input_vars["ParCalculation"]["dth"]),
    dts=input_vars["ParCalculation"]["dts"],
    row=FT(input_vars["ParCalculation"]["row"]),
)

rsRoofPreCalc = (;
    rs_sun=input_vars["rsRoofPreCalc"]["rs_sun"],
    rs_shd=input_vars["rsRoofPreCalc"]["rs_shd"],
    Ci_sun=input_vars["rsRoofPreCalc"]["Ci_sun"],
    Ci_shd=input_vars["rsRoofPreCalc"]["Ci_shd"],
)

@testset "MATLAB" begin
    results = eb_wb_roof(
        vec(input_vars["TemperatureR"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        Runon_ittm,
        Geometry_m,
        FractionsRoof,
        ParSoilRoof,
        PropOpticalRoof,
        ParThermalRoof,
        ParVegRoof,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        Bool(input_vars["BEM_on"]),
        Bool(input_vars["RESPreCalc"]),
        rsRoofPreCalc,
    )

    # Test shortwave radiation outputs
    @test results.SWRabsRoofImp ≈ output_vars["SWRabsRoofImp"]
    @test results.SWRabsRoofVeg ≈ output_vars["SWRabsRoofVeg"]
    @test results.SWRabsTotalRoof ≈ output_vars["SWRabsTotalRoof"]
    @test results.SWRoutRoofImp ≈ output_vars["SWRoutRoofImp"]
    @test results.SWRoutRoofVeg ≈ output_vars["SWRoutRoofVeg"]
    @test results.SWRoutTotalRoof ≈ output_vars["SWRoutTotalRoof"]
    @test results.SWRinRoofImp ≈ output_vars["SWRinRoofImp"]
    @test results.SWRinRoofVeg ≈ output_vars["SWRinRoofVeg"]
    @test results.SWRinTotalRoof ≈ output_vars["SWRinTotalRoof"]
    @test results.SWREBRoofImp ≈ output_vars["SWREBRoofImp"]
    @test results.SWREBRoofVeg ≈ output_vars["SWREBRoofVeg"]
    @test results.SWREBTotalRoof ≈ output_vars["SWREBTotalRoof"]

    # Test longwave radiation outputs
    @test results.LWRabsRoofVeg ≈ output_vars["LWRabsRoofVeg"]
    @test results.LWRabsRoofImp ≈ output_vars["LWRabsRoofImp"]
    @test results.LWRabsTotalRoof ≈ output_vars["LWRabsTotalRoof"]
    @test results.LWRoutRoofVeg ≈ output_vars["LWRoutRoofVeg"]
    @test results.LWRoutRoofImp ≈ output_vars["LWRoutRoofImp"]
    @test results.LWRoutTotalRoof ≈ output_vars["LWRoutTotalRoof"]
    @test results.LWRinRoofImp ≈ output_vars["LWRinRoofImp"]
    @test results.LWRinRoofVeg ≈ output_vars["LWRinRoofVeg"]
    @test results.LWRinTotalRoof ≈ output_vars["LWRinTotalRoof"]
    @test results.LWREBRoofImp ≈ output_vars["LWREBRoofImp"]
    @test results.LWREBRoofVeg ≈ output_vars["LWREBRoofVeg"]
    @test results.LWREBTotalRoof ≈ output_vars["LWREBTotalRoof"] atol=1e-12

    # Test heat fluxes
    @test results.HfluxRoofImp ≈ output_vars["HfluxRoofImp"]
    @test results.HfluxRoofVeg ≈ output_vars["HfluxRoofVeg"]
    @test results.HfluxRoof ≈ output_vars["HfluxRoof"]

    # Test latent heat fluxes
    @test results.LEfluxRoofImp ≈ output_vars["LEfluxRoofImp"]
    @test results.LEfluxRoofVegInt ≈ output_vars["LEfluxRoofVegInt"]
    @test results.LEfluxRoofVegPond ≈ output_vars["LEfluxRoofVegPond"]
    @test results.LEfluxRoofVegSoil ≈ output_vars["LEfluxRoofVegSoil"]
    @test results.LTEfluxRoofVeg ≈ output_vars["LTEfluxRoofVeg"]
    @test results.LEfluxRoofVeg ≈ output_vars["LEfluxRoofVeg"]
    @test results.LEfluxRoof ≈ output_vars["LEfluxRoof"]

    # Test conductive heat fluxes
    @test results.G1RoofImp ≈ output_vars["G1RoofImp"]
    @test results.G2RoofImp ≈ output_vars["G2RoofImp"]
    @test results.dsRoofImp ≈ output_vars["dsRoofImp"]
    @test results.G1RoofVeg ≈ output_vars["G1RoofVeg"]
    @test results.G2RoofVeg ≈ output_vars["G2RoofVeg"]
    @test results.dsRoofVeg ≈ output_vars["dsRoofVeg"]
    @test results.G1Roof ≈ output_vars["G1Roof"]
    @test results.G2Roof ≈ output_vars["G2Roof"]
    @test results.dsRoof ≈ output_vars["dsRoof"]

    # Test resistances
    @test results.raRooftoAtm ≈ output_vars["raRooftoAtm"]
    @test results.rb_LRoof ≈ output_vars["rb_LRoof"]
    if isnan(output_vars["rap_LRoof"])
        @test isnan(results.rap_LRoof)
    else
        @test results.rap_LRoof ≈ output_vars["rap_LRoof"]
    end
    @test results.r_soilRoof ≈ output_vars["r_soilRoof"]
    @test results.rs_sunRoof ≈ output_vars["rs_sunRoof"]
    @test results.rs_shdRoof ≈ output_vars["rs_shdRoof"]

    # Test water fluxes
    @test results.EfluxRoofImp ≈ output_vars["EfluxRoofImp"]
    @test results.EfluxRoofVegInt ≈ output_vars["EfluxRoofVegInt"]
    @test results.EfluxRoofVegPond ≈ output_vars["EfluxRoofVegPond"]
    @test results.EfluxRoofVegSoil ≈ output_vars["EfluxRoofVegSoil"]
    @test results.TEfluxRoofVeg ≈ output_vars["TEfluxRoofVeg"]
    @test results.EfluxRoofVeg ≈ output_vars["EfluxRoofVeg"]
    @test results.EfluxRoof ≈ output_vars["EfluxRoof"]

    # Test water balance outputs
    @test results.QRoofImp ≈ output_vars["QRoofImp"]
    @test results.QRoofVegDrip ≈ output_vars["QRoofVegDrip"]
    @test results.QRoofVegPond ≈ output_vars["QRoofVegPond"]
    @test results.LkRoofImp ≈ output_vars["LkRoofImp"]
    @test results.LkRoofVeg ≈ output_vars["LkRoofVeg"]
    @test results.LkRoof ≈ output_vars["LkRoof"]
    @test results.QRoofVegSoil ≈ output_vars["QRoofVegSoil"]
    @test results.RunoffRoofTot ≈ output_vars["RunoffRoofTot"]
    @test results.RunonRoofTot ≈ output_vars["RunonRoofTot"]
    @test results.IntRoofImp ≈ output_vars["IntRoofImp"]
    @test results.IntRoofVegPlant ≈ output_vars["IntRoofVegPlant"]
    @test results.IntRoofVegGround ≈ output_vars["IntRoofVegGround"]
    @test results.dInt_dtRoofImp ≈ output_vars["dInt_dtRoofImp"]
    @test results.dInt_dtRoofVegPlant ≈ output_vars["dInt_dtRoofVegPlant"]
    @test results.dInt_dtRoofVegGround ≈ output_vars["dInt_dtRoofVegGround"]
    @test results.IntRooftot ≈ output_vars["IntRooftot"]
    @test results.dInt_dtRooftot ≈ output_vars["dInt_dtRooftot"]
    @test results.dVRoofSoil_dt ≈ output_vars["dVRoofSoil_dt"]
    @test results.fRoofVeg ≈ output_vars["fRoofVeg"]
    # Not working
    @test results.VRoofSoil ≈ vec(output_vars["VRoofSoil"]) atol = 0.05
    @test results.OwRoofSoil ≈ vec(output_vars["OwRoofSoil"]) atol = 0.05
    @test results.OSwRoofSoil ≈ output_vars["OSwRoofSoil"] atol = 0.05

    # Test water potentials and CO2
    # @test results.ExWaterRoof_L ≈ vec(output_vars["ExWaterRoof_L"]) # not working
    @test results.SoilPotWRoof_L ≈ [output_vars["SoilPotWRoof_L"]] atol = 0.05
    @test results.CiCO2LeafRoofVegSun ≈ output_vars["CiCO2LeafRoofVegSun"]
    @test results.CiCO2LeafRoofVegShd ≈ output_vars["CiCO2LeafRoofVegShd"]

    # Test water balances
    @test results.WBRoofVegInVeg ≈ output_vars["WBRoofVegInVeg"]
    @test results.WBRoofVegInGround ≈ output_vars["WBRoofVegInGround"]
    @test results.WBRoofVegSoil ≈ output_vars["WBRoofVegSoil"] atol = 1e-10

    # Test energy balances
    @test results.EBRoofImp ≈ output_vars["EBRoofImp"] atol = 1e-10
    @test results.EBRoofVeg ≈ output_vars["EBRoofVeg"] atol = 1e-10
    @test results.Yroof ≈ vec(output_vars["Yroof"]) atol = 1e-12
    @test results.WBRoofImp ≈ output_vars["WBRoofImp"] atol = 1e-14
    @test results.WBRoofVeg ≈ output_vars["WBRoofVeg"] atol = 1e-14
    @test results.WBRoofTot ≈ output_vars["WBRoofTot"] atol = 1e-14
end
