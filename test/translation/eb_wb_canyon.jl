using Test
using MAT
using UrbanTethysChloris: eb_wb_canyon
using UrbanTethysChloris.RayTracing: ViewFactor
using UrbanTethysChloris.ModelComponents.Parameters
using UrbanTethysChloris.ModelComponents.ForcingInputs
using UrbanTethysChloris.ModelComponents.ModelVariables
using ...TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("EB_WB_canyon.json")

view_factor = ViewFactor(FT, input_vars["ViewFactor"])
Gemeotry_m = UrbanGeometryParameters(
    FT, input_vars["Gemeotry_m"], input_vars["geometry"], input_vars["ParTree"]
)
FractionsGround = LocationSpecificSurfaceFractions(FT, input_vars["FractionsGround"])
ParSoilGround = VegetatedSoilParameters(FT, input_vars["ParSoilGround"], "_G")
PropOpticalGround = VegetatedOpticalProperties(FT, input_vars["PropOpticalGround"])
PropOpticalWall = SimpleOpticalProperties(FT, input_vars["PropOpticalWall"])
PropOpticalTree = SimpleOpticalProperties(FT, input_vars["PropOpticalTree"])
ParThermalGround = LocationSpecificThermalProperties(
    FT, input_vars["ParThermalGround"], "_imp"
)
ParThermalWall = LocationSpecificThermalProperties(FT, input_vars["ParThermalWall"])
ParVegGround = HeightDependentVegetationParameters(FT, input_vars["ParVegGround"])
ParVegTree = HeightDependentVegetationParameters(FT, input_vars["ParVegTree"])
PropOpticalIndoors = IndoorOpticalProperties(FT, input_vars["PropOpticalIndoors"])
ParHVAC = HVACParameters(FT, input_vars["ParHVAC"])
ParWindows = WindowParameters(FT, input_vars["ParWindows"])
ParThermalBulidFloor = ThermalBuilding(FT, input_vars["ParThermalBulidFloor"])

ParInterceptionTree = (; Sp_In=input_vars["ParInterceptionTree"]["Sp_In"],)
@testset "Model variables" begin
    TempVec_ittm = TempVec(FT, input_vars["TempVec_ittm"])
    Humidity_ittm = Humidity(FT, input_vars["Humidity_ittm"])
    MeteoData = MeteorologicalInputs(
        FT,
        input_vars["MeteoData"],
        input_vars["HumidityAtm"],
        input_vars["ParCalculation"]["cp_atm"],
        input_vars["ParCalculation"]["rho_atm"],
    )
    Int_ittm = Interception(FT, input_vars["Int_ittm"])
    ExWater_ittm = ExWater(FT, input_vars["ExWater_ittm"])
    Vwater_ittm = Vwater(FT, input_vars["Vwater_ittm"])
    Owater_ittm = Owater(FT, input_vars["Owater_ittm"])
    SoilPotW_ittm = SoilPotW(FT, input_vars["SoilPotW_ittm"])
    CiCO2Leaf_ittm = CiCO2Leaf(FT, input_vars["CiCO2Leaf_ittm"])
    TempDamp_ittm = TempDamp(FT, input_vars["TempDamp_ittm"])
    Runon_ittm = Runon(FT, input_vars["Runon_ittm"])
    Qinlat_ittm = Qinlat(FT, input_vars["Qinlat_ittm"])
    TempVecB_ittm = TempVecB(FT, input_vars["TempVecB_ittm"])

    HumidityAtm = MeteorologicalInputs(
        FT,
        input_vars["MeteoData"],
        input_vars["HumidityAtm"],
        input_vars["ParCalculation"]["cp_atm"],
        input_vars["ParCalculation"]["rho_atm"],
    )

    Anthropogenic = AnthropogenicInputs(FT, input_vars["Anthropogenic"])

    ParCalculation = (;
        dth=FT(input_vars["ParCalculation"]["dth"]),
        dts=Int(input_vars["ParCalculation"]["dts"]),
        row=FT(input_vars["ParCalculation"]["row"]),
    )

    HVAC_Schedule = HVACSchedule(FT, input_vars["HVACSchedule"])

    fconvPreCalc = input_vars["fconvPreCalc"]

    fconv = input_vars["fconv"]

    rsGroundPreCalc = (;
        rs_sun_L=input_vars["rsGroundPreCalc"]["rs_sun_L"],
        rs_shd_L=input_vars["rsGroundPreCalc"]["rs_shd_L"],
        Ci_sun_L=input_vars["rsGroundPreCalc"]["Ci_sun_L"],
        Ci_shd_L=input_vars["rsGroundPreCalc"]["Ci_shd_L"],
    )

    rsTreePreCalc = (;
        rs_sun_H=input_vars["rsTreePreCalc"]["rs_sun_H"],
        rs_shd_H=input_vars["rsTreePreCalc"]["rs_shd_H"],
        Ci_sun_H=input_vars["rsTreePreCalc"]["Ci_sun_H"],
        Ci_shd_H=input_vars["rsTreePreCalc"]["Ci_shd_H"],
    )

    WallLayers = (;
        dz1_wall=input_vars["WallLayers"]["dz1_wall"],
        dz2_wall=input_vars["WallLayers"]["dz2_wall"],
    )
    Sun_Position = SunPositionInputs(FT, input_vars["SunPosition"])

    results = eb_wb_canyon(
        vec(input_vars["TemperatureC"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
        Humidity_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        TempDamp_ittm,
        Runon_ittm,
        Qinlat_ittm,
        view_factor,
        Gemeotry_m,
        FractionsGround,
        WallLayers,
        ParSoilGround,
        ParInterceptionTree,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParThermalGround,
        ParThermalWall,
        ParVegGround,
        ParVegTree,
        Sun_Position,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        TempVecB_ittm,
        FT(input_vars["G2Roof"]),
        PropOpticalIndoors,
        ParHVAC,
        ParThermalBulidFloor,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        Bool(input_vars["RESPreCalc"]),
        Bool(fconvPreCalc),
        FT(fconv),
        rsGroundPreCalc,
        rsTreePreCalc,
        HVAC_Schedule,
    )

    @test results.Ycanyon ≈ vec(output_vars["Ycanyon"]) atol=1e-6
    @test results.G2WallSun ≈ output_vars["G2WallSun"]
    @test results.G2WallShade ≈ output_vars["G2WallShade"]

    @test results.SWRabs_t.GroundImp ≈ output_vars["SWRabs_t"]["SWRabsGroundImp"]
    @test results.SWRabs_t.GroundBare ≈ output_vars["SWRabs_t"]["SWRabsGroundBare"]
    @test results.SWRabs_t.GroundVeg ≈ output_vars["SWRabs_t"]["SWRabsGroundVeg"]
    @test results.SWRabs_t.WallSun ≈ output_vars["SWRabs_t"]["SWRabsWallSun"]
    @test results.SWRabs_t.WallShade ≈ output_vars["SWRabs_t"]["SWRabsWallShade"]
    @test results.SWRabs_t.Tree ≈ output_vars["SWRabs_t"]["SWRabsTree"]
    @test results.SWRabs_t.TotalGround ≈ output_vars["SWRabs_t"]["SWRabsTotalGround"]
    @test results.SWRabs_t.TotalCanyon ≈ output_vars["SWRabs_t"]["SWRabsTotalCanyon"]
end
@testset "Named tuples" begin

    # Create named tuples for other parameters
    TempVec_ittm = (;
        TWallSun=input_vars["TempVec_ittm"]["TWallSun"],
        TWallShade=input_vars["TempVec_ittm"]["TWallShade"],
        TWallIntSun=input_vars["TempVec_ittm"]["TWallIntSun"],
        TWallIntShade=input_vars["TempVec_ittm"]["TWallIntShade"],
        TGroundImp=input_vars["TempVec_ittm"]["TGroundImp"],
        TGroundBare=input_vars["TempVec_ittm"]["TGroundBare"],
        TGroundVeg=input_vars["TempVec_ittm"]["TGroundVeg"],
        TTree=input_vars["TempVec_ittm"]["TTree"],
        TCanyon=input_vars["TempVec_ittm"]["TCanyon"],
        T2m=input_vars["TempVec_ittm"]["T2m"],
    )

    Humidity_ittm = (;
        CanyonSpecific=input_vars["Humidity_ittm"]["CanyonSpecific"],
        q2m=input_vars["Humidity_ittm"]["q2m"],
    )

    MeteoData = (;
        SW_dir=FT(input_vars["MeteoData"]["SW_dir"]),
        SW_diff=FT(input_vars["MeteoData"]["SW_diff"]),
        LWR=input_vars["MeteoData"]["LWR"],
        Rain=FT(input_vars["MeteoData"]["Rain"]),
        Tatm=input_vars["MeteoData"]["Tatm"],
        Pre=input_vars["MeteoData"]["Pre"],
        ea=input_vars["MeteoData"]["ea"],
        Zatm=FT(input_vars["MeteoData"]["Zatm"]),
        Uatm=input_vars["MeteoData"]["Uatm"],
        q_atm=input_vars["MeteoData"]["q_atm"],
        Catm_O2=FT(input_vars["MeteoData"]["Catm_O2"]),
        Catm_CO2=FT(input_vars["MeteoData"]["Catm_CO2"]),
    )

    Int_ittm = (;
        IntGroundImp=FT(input_vars["Int_ittm"]["IntGroundImp"]),
        IntGroundVegPlant=FT(input_vars["Int_ittm"]["IntGroundVegPlant"]),
        IntGroundVegGround=FT(input_vars["Int_ittm"]["IntGroundVegGround"]),
        IntTree=FT(input_vars["Int_ittm"]["IntTree"]),
        IntGroundBare=FT(input_vars["Int_ittm"]["IntGroundBare"]),
    )

    ExWater_ittm = (;
        ExWaterGroundImp_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundImp_H"]),
        ExWaterGroundBare_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundBare_H"]),
        ExWaterGroundVeg_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_H"]),
        ExWaterGroundVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_L"]),
    )

    Vwater_ittm = (;
        VGroundSoilVeg=vec(input_vars["Vwater_ittm"]["VGroundSoilVeg"]),
        VGroundSoilImp=vec(input_vars["Vwater_ittm"]["VGroundSoilImp"]),
        VGroundSoilBare=vec(input_vars["Vwater_ittm"]["VGroundSoilBare"]),
    )

    Owater_ittm = (;
        OwGroundSoilVeg=vec(input_vars["Owater_ittm"]["OwGroundSoilVeg"]),
        OwGroundSoilBare=vec(input_vars["Owater_ittm"]["OwGroundSoilBare"]),
        OwGroundSoilImp=vec(input_vars["Owater_ittm"]["OwGroundSoilImp"]),
    )

    SoilPotW_ittm = (;
        SoilPotWGroundVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWGroundVeg_L"],
        SoilPotWGroundTot_H=input_vars["SoilPotW_ittm"]["SoilPotWGroundTot_H"],
    )

    CiCO2Leaf_ittm = (;
        CiCO2LeafGroundVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegSun"],
        CiCO2LeafGroundVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegShd"],
        CiCO2LeafTreeSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeSun"],
        CiCO2LeafTreeShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeShd"],
    )

    TempDamp_ittm = (;
        TDampGroundImp=input_vars["TempDamp_ittm"]["TDampGroundImp"],
        TDampGroundBare=input_vars["TempDamp_ittm"]["TDampGroundBare"],
        TDampGroundVeg=input_vars["TempDamp_ittm"]["TDampGroundVeg"],
        TDampGroundBuild=input_vars["TempDamp_ittm"]["TDampGroundBuild"],
    )

    Runon_ittm = (; RunonGroundTot=FT(input_vars["Runon_ittm"]["RunonGroundTot"]),)

    Qinlat_ittm = (;
        Qin_imp=FT.(vec(input_vars["Qinlat_ittm"]["Qin_imp"])),
        Qin_bare=FT.(vec(input_vars["Qinlat_ittm"]["Qin_bare"])),
        Qin_veg=FT.(vec(input_vars["Qinlat_ittm"]["Qin_veg"])),
    )

    TempVecB_ittm = (;
        Tbin=input_vars["TempVecB_ittm"]["Tbin"],
        qbin=input_vars["TempVecB_ittm"]["qbin"],
        Tinground=input_vars["TempVecB_ittm"]["Tinground"],
        Tintmass=input_vars["TempVecB_ittm"]["Tintmass"],
        Tinwallsun=input_vars["TempVecB_ittm"]["Tinwallsun"],
        Tinwallshd=input_vars["TempVecB_ittm"]["Tinwallshd"],
        Tceiling=input_vars["TempVecB_ittm"]["Tceiling"],
    )

    HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])

    Anthropogenic = (;
        Qf_canyon=input_vars["Anthropogenic"]["Qf_canyon"],
        Waterf_canyonBare=FT(input_vars["Anthropogenic"]["Waterf_canyonBare"]),
        Waterf_canyonVeg=FT(input_vars["Anthropogenic"]["Waterf_canyonVeg"]),
    )

    ParCalculation = (;
        dth=FT(input_vars["ParCalculation"]["dth"]),
        dts=Int(input_vars["ParCalculation"]["dts"]),
        row=FT(input_vars["ParCalculation"]["row"]),
    )

    HVAC_Schedule = (;
        Hequip=FT(input_vars["HVACSchedule"]["Hequip"]),
        Hpeople=FT(input_vars["HVACSchedule"]["Hpeople"]),
        LEequip=FT(input_vars["HVACSchedule"]["LEequip"]),
        LEpeople=FT(input_vars["HVACSchedule"]["LEpeople"]),
        AirConRoomFraction=FT(input_vars["HVACSchedule"]["AirConRoomFraction"]),
    )

    fconvPreCalc = input_vars["fconvPreCalc"]

    fconv = input_vars["fconv"]

    rsGroundPreCalc = (;
        rs_sun_L=input_vars["rsGroundPreCalc"]["rs_sun_L"],
        rs_shd_L=input_vars["rsGroundPreCalc"]["rs_shd_L"],
        Ci_sun_L=input_vars["rsGroundPreCalc"]["Ci_sun_L"],
        Ci_shd_L=input_vars["rsGroundPreCalc"]["Ci_shd_L"],
    )

    rsTreePreCalc = (;
        rs_sun_H=input_vars["rsTreePreCalc"]["rs_sun_H"],
        rs_shd_H=input_vars["rsTreePreCalc"]["rs_shd_H"],
        Ci_sun_H=input_vars["rsTreePreCalc"]["Ci_sun_H"],
        Ci_shd_H=input_vars["rsTreePreCalc"]["Ci_shd_H"],
    )

    WallLayers = (;
        dz1_wall=input_vars["WallLayers"]["dz1_wall"],
        dz2_wall=input_vars["WallLayers"]["dz2_wall"],
    )

    Sun_Position = (;
        theta_n=input_vars["SunPosition"]["theta_n"],
        theta_Z=input_vars["SunPosition"]["theta_Z"],
    )

    results = eb_wb_canyon(
        vec(input_vars["TemperatureC"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
        Humidity_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        TempDamp_ittm,
        Runon_ittm,
        Qinlat_ittm,
        view_factor,
        Gemeotry_m,
        FractionsGround,
        WallLayers,
        ParSoilGround,
        ParInterceptionTree,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParThermalGround,
        ParThermalWall,
        ParVegGround,
        ParVegTree,
        Sun_Position,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        TempVecB_ittm,
        FT(input_vars["G2Roof"]),
        PropOpticalIndoors,
        ParHVAC,
        ParThermalBulidFloor,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        Bool(input_vars["RESPreCalc"]),
        FT(fconvPreCalc),
        FT(fconv),
        rsGroundPreCalc,
        rsTreePreCalc,
        HVAC_Schedule,
    )

    @test results.Ycanyon ≈ vec(output_vars["Ycanyon"]) atol=1e-6
    @test results.G2WallSun ≈ output_vars["G2WallSun"]
    @test results.G2WallShade ≈ output_vars["G2WallShade"]

    @test results.SWRabs_t.GroundImp ≈ output_vars["SWRabs_t"]["SWRabsGroundImp"]
    @test results.SWRabs_t.GroundBare ≈ output_vars["SWRabs_t"]["SWRabsGroundBare"]
    @test results.SWRabs_t.GroundVeg ≈ output_vars["SWRabs_t"]["SWRabsGroundVeg"]
    @test results.SWRabs_t.WallSun ≈ output_vars["SWRabs_t"]["SWRabsWallSun"]
    @test results.SWRabs_t.WallShade ≈ output_vars["SWRabs_t"]["SWRabsWallShade"]
    @test results.SWRabs_t.Tree ≈ output_vars["SWRabs_t"]["SWRabsTree"]
    @test results.SWRabs_t.TotalGround ≈ output_vars["SWRabs_t"]["SWRabsTotalGround"]
    @test results.SWRabs_t.TotalCanyon ≈ output_vars["SWRabs_t"]["SWRabsTotalCanyon"]
end
