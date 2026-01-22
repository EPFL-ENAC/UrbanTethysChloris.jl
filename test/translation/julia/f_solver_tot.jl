using Test
using UrbanTethysChloris: f_solver_tot
using UrbanTethysChloris: ExtrapolatedTempVec, ExtrapolatedHumidity, ExtrapolatedTempVecB
using UrbanTethysChloris: Meteotm1
using UrbanTethysChloris.RayTracing: ViewFactor
using ...TestUtils: load_matlab_data
using UrbanTethysChloris.ModelComponents.Parameters
using UrbanTethysChloris.ModelComponents.ModelVariables
using UrbanTethysChloris.ModelComponents.ForcingInputs

FT = Float64
input_vars, output_vars = load_matlab_data("fSolver_Tot.json")

view_factor = ViewFactor(FT, input_vars["ViewFactor"])
Gemeotry_m = UrbanGeometryParameters(
    FT, input_vars["Gemeotry_m"], input_vars["geometry"], input_vars["ParTree"]
)
FractionsGround = LocationSpecificSurfaceFractions(FT, input_vars["FractionsGround"])
FractionsRoof = LocationSpecificSurfaceFractions(FT, input_vars["FractionsRoof"])
ParSoilGround = VegetatedSoilParameters(FT, input_vars["ParSoilGround"], "_G")
ParSoilRoof = VegetatedSoilParameters(FT, input_vars["ParSoilRoof"], "_R")
PropOpticalGround = VegetatedOpticalProperties(FT, input_vars["PropOpticalGround"])
PropOpticalRoof = VegetatedOpticalProperties(FT, input_vars["PropOpticalRoof"])
PropOpticalWall = SimpleOpticalProperties(FT, input_vars["PropOpticalWall"])
PropOpticalTree = SimpleOpticalProperties(FT, input_vars["PropOpticalTree"])
ParThermalGround = LocationSpecificThermalProperties(
    FT, input_vars["ParThermalGround"], "_imp"
)
ParThermalWall = LocationSpecificThermalProperties(FT, input_vars["ParThermalWall"])
ParThermalRoof = LocationSpecificThermalProperties(FT, input_vars["ParThermalRoof"], "_imp")
ParVegGround = HeightDependentVegetationParameters(FT, input_vars["ParVegGround"])
ParVegTree = HeightDependentVegetationParameters(FT, input_vars["ParVegTree"])
ParVegRoof = HeightDependentVegetationParameters(FT, input_vars["ParVegRoof"])
PropOpticalIndoors = IndoorOpticalProperties(FT, input_vars["PropOpticalIndoors"])
ParHVAC = HVACParameters(FT, input_vars["ParHVAC"])
ParThermalBulidFloor = ThermalBuilding(FT, input_vars["ParThermalBulidFloor"])
ParWindows = WindowParameters(FT, input_vars["ParWindows"])

ParInterceptionTree = (; Sp_In=input_vars["ParInterceptionTree"]["Sp_In"],)
WallLayers = (;
    dz1_wall=input_vars["WallLayers"]["dz1_wall"],
    dz2_wall=input_vars["WallLayers"]["dz2_wall"],
)
ParCalculation = (;
    dth=Int(input_vars["ParCalculation"]["dth"]),
    dts=Int(input_vars["ParCalculation"]["dts"]),
    row=FT(input_vars["ParCalculation"]["row"]),
)

rsRoofPreCalc = (;
    rs_sun=input_vars["rsRoofPreCalc"]["rs_sun"],
    rs_shd=input_vars["rsRoofPreCalc"]["rs_shd"],
    Ci_sun=input_vars["rsRoofPreCalc"]["Ci_sun"],
    Ci_shd=input_vars["rsRoofPreCalc"]["Ci_shd"],
)

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

@testset "Composite types" begin
    TempVec_ittm = TempVec(FT, input_vars["TempVec_ittm"])
    TempVecB_ittm = TempVecB(FT, input_vars["TempVecB_ittm"])
    Humidity_ittm = Humidity(FT, input_vars["Humidity_ittm"])
    Int_ittm = Interception(FT, input_vars["Int_ittm"])
    ExWater_ittm = ExWater(FT, input_vars["ExWater_ittm"])
    Vwater_ittm = Vwater(FT, input_vars["Vwater_ittm"])
    Owater_ittm = Owater(FT, input_vars["Owater_ittm"])
    SoilPotW_ittm = SoilPotW(FT, input_vars["SoilPotW_ittm"])
    CiCO2Leaf_ittm = CiCO2Leaf(FT, input_vars["CiCO2Leaf_ittm"])
    TempDamp_ittm = TempDamp(FT, input_vars["TempDamp_ittm"])
    MeteoData = MeteorologicalInputs(
        FT,
        input_vars["MeteoData"],
        input_vars["HumidityAtm"],
        input_vars["ParCalculation"]["cp_atm"],
        input_vars["ParCalculation"]["rho_atm"],
    )
    HumidityAtm = MeteorologicalInputs(
        FT,
        input_vars["MeteoData"],
        input_vars["HumidityAtm"],
        input_vars["ParCalculation"]["cp_atm"],
        input_vars["ParCalculation"]["rho_atm"],
    )
    SunPosition = SunPositionInputs(FT, input_vars["SunPosition"])
    Anthropogenic = AnthropogenicInputs(FT, input_vars["Anthropogenic"])
    HVAC_Schedule = HVACSchedule(FT, input_vars["HVACSchedule"])

    TempVec_ittm2Ext = ExtrapolatedTempVec{FT}(
        TWallSun=input_vars["TempVec_ittm2Ext"]["TWallSun"],
        TWallShade=input_vars["TempVec_ittm2Ext"]["TWallShade"],
        TWallIntSun=input_vars["TempVec_ittm2Ext"]["TWallIntSun"],
        TWallIntShade=input_vars["TempVec_ittm2Ext"]["TWallIntShade"],
        TGroundImp=input_vars["TempVec_ittm2Ext"]["TGroundImp"],
        TGroundBare=input_vars["TempVec_ittm2Ext"]["TGroundBare"],
        TGroundVeg=input_vars["TempVec_ittm2Ext"]["TGroundVeg"],
        TTree=input_vars["TempVec_ittm2Ext"]["TTree"],
        TCanyon=input_vars["TempVec_ittm2Ext"]["TCanyon"],
        TRoofVeg=input_vars["TempVec_ittm2Ext"]["TRoofVeg"],
        TRoofIntVeg=input_vars["TempVec_ittm2Ext"]["TRoofIntVeg"],
        TRoofIntImp=input_vars["TempVec_ittm2Ext"]["TRoofIntImp"],
        TRoofImp=input_vars["TempVec_ittm2Ext"]["TRoofImp"],
        Tatm=fill(zero(FT), 4),
        T2m=fill(zero(FT), 4),
    )
    Humidity_ittm2Ext = ExtrapolatedHumidity{FT}(
        CanyonSpecific=input_vars["Humidity_ittm2Ext"]["CanyonSpecific"],
        CanyonRelative=fill(zero(FT), 4),
        CanyonVapourPre=fill(zero(FT), 4),
        CanyonRelativeSat=fill(zero(FT), 4),
        CanyonSpecificSat=fill(zero(FT), 4),
        CanyonVapourPreSat=fill(zero(FT), 4),
        AtmRelative=fill(zero(FT), 4),
        AtmSpecific=fill(zero(FT), 4),
        AtmVapourPre=fill(zero(FT), 4),
        AtmRelativeSat=fill(zero(FT), 4),
        AtmSpecificSat=fill(zero(FT), 4),
        AtmVapourPreSat=fill(zero(FT), 4),
        q2m=fill(zero(FT), 4),
    )
    TempVecB_ittm2Ext = ExtrapolatedTempVecB{FT}(
        Tbin=input_vars["TempVecB_ittm2Ext"]["Tbin"],
        qbin=input_vars["TempVecB_ittm2Ext"]["qbin"],
        Tinground=input_vars["TempVecB_ittm2Ext"]["Tinground"],
        Tintmass=input_vars["TempVecB_ittm2Ext"]["Tintmass"],
        Tinwallsun=input_vars["TempVecB_ittm2Ext"]["Tinwallsun"],
        Tinwallshd=input_vars["TempVecB_ittm2Ext"]["Tinwallshd"],
        Tceiling=input_vars["TempVecB_ittm2Ext"]["Tceiling"],
        Twindows=input_vars["TempVecB_ittm2Ext"]["Twindows"],
    )
    Meteo_ittm = Meteotm1{FT}(
        SWRin=FT.(input_vars["Meteo_ittm"]["SWRin"]),
        Rain=FT.(input_vars["Meteo_ittm"]["Rain"]),
    )

    T, fval, exitflag = f_solver_tot(
        TempVec_ittm,
        TempVecB_ittm,
        Humidity_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        TempDamp_ittm,
        view_factor,
        Gemeotry_m,
        FractionsGround,
        FractionsRoof,
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
        ParSoilRoof,
        PropOpticalRoof,
        ParThermalRoof,
        ParVegRoof,
        SunPosition,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        PropOpticalIndoors,
        ParHVAC,
        ParThermalBulidFloor,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        TempVec_ittm2Ext,
        Humidity_ittm2Ext,
        TempVecB_ittm2Ext,
        Meteo_ittm,
        Bool(input_vars["RESPreCalc"]),
        Bool(input_vars["fconvPreCalc"]),
        FT(input_vars["fconv"]),
        rsRoofPreCalc,
        rsGroundPreCalc,
        rsTreePreCalc,
        HVAC_Schedule,
    )

    @test T ≈ vec(output_vars["T"]) atol=1e-10
    @test fval ≈ vec(output_vars["fval"]) atol=1e-10
    @test exitflag == output_vars["exitflag"]
end

@testset "Named tuples" begin
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
        TRoofVeg=input_vars["TempVec_ittm"]["TRoofVeg"],
        TRoofIntVeg=input_vars["TempVec_ittm"]["TRoofIntVeg"],
        TRoofIntImp=input_vars["TempVec_ittm"]["TRoofIntImp"],
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

    Humidity_ittm = (; CanyonSpecific=input_vars["Humidity_ittm"]["CanyonSpecific"])

    MeteoData = (;
        SW_dir=FT(input_vars["MeteoData"]["SW_dir"]),
        SW_diff=FT(input_vars["MeteoData"]["SW_diff"]),
        LWR=input_vars["MeteoData"]["LWR"],
        Rain=input_vars["MeteoData"]["Rain"],
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
        IntGroundImp=input_vars["Int_ittm"]["IntGroundImp"],
        IntGroundVegPlant=input_vars["Int_ittm"]["IntGroundVegPlant"],
        IntGroundVegGround=input_vars["Int_ittm"]["IntGroundVegGround"],
        IntTree=input_vars["Int_ittm"]["IntTree"],
        IntGroundBare=input_vars["Int_ittm"]["IntGroundBare"],
        IntRoofImp=input_vars["Int_ittm"]["IntRoofImp"],
        IntRoofVegPlant=input_vars["Int_ittm"]["IntRoofVegPlant"],
        IntRoofVegGround=input_vars["Int_ittm"]["IntRoofVegGround"],
    )

    ExWater_ittm = (;
        ExWaterGroundImp_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundImp_H"]),
        ExWaterGroundBare_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundBare_H"]),
        ExWaterGroundVeg_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_H"]),
        ExWaterGroundVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_L"]),
        ExWaterRoofVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterRoofVeg_L"]),
    )

    Vwater_ittm = (;
        VGroundSoilVeg=vec(input_vars["Vwater_ittm"]["VGroundSoilVeg"]),
        VGroundSoilImp=vec(input_vars["Vwater_ittm"]["VGroundSoilImp"]),
        VGroundSoilBare=vec(input_vars["Vwater_ittm"]["VGroundSoilBare"]),
        VRoofSoilVeg=vec(input_vars["Vwater_ittm"]["VRoofSoilVeg"]),
    )

    Owater_ittm = (;
        OwGroundSoilVeg=vec(input_vars["Owater_ittm"]["OwGroundSoilVeg"]),
        OwGroundSoilBare=vec(input_vars["Owater_ittm"]["OwGroundSoilBare"]),
        OwGroundSoilImp=vec(input_vars["Owater_ittm"]["OwGroundSoilImp"]),
        OwRoofSoilVeg=vec(input_vars["Owater_ittm"]["OwRoofSoilVeg"]),
    )

    SoilPotW_ittm = (;
        SoilPotWGroundVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWGroundVeg_L"],
        SoilPotWGroundTot_H=input_vars["SoilPotW_ittm"]["SoilPotWGroundTot_H"],
        SoilPotWRoofVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWRoofVeg_L"],
    )

    CiCO2Leaf_ittm = (;
        CiCO2LeafGroundVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegSun"],
        CiCO2LeafGroundVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegShd"],
        CiCO2LeafTreeSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeSun"],
        CiCO2LeafTreeShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeShd"],
        CiCO2LeafRoofVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegSun"],
        CiCO2LeafRoofVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegShd"],
    )

    TempDamp_ittm = (;
        TDampGroundImp=input_vars["TempDamp_ittm"]["TDampGroundImp"],
        TDampGroundBare=input_vars["TempDamp_ittm"]["TDampGroundBare"],
        TDampGroundVeg=input_vars["TempDamp_ittm"]["TDampGroundVeg"],
        TDampGroundBuild=input_vars["TempDamp_ittm"]["TDampGroundBuild"],
    )

    SunPosition = (;
        theta_n=input_vars["SunPosition"]["theta_n"],
        theta_Z=input_vars["SunPosition"]["theta_Z"],
    )

    HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])

    Anthropogenic = (;
        Qf_canyon=input_vars["Anthropogenic"]["Qf_canyon"],
        Tb=input_vars["Anthropogenic"]["Tb"],
        Waterf_roof=input_vars["Anthropogenic"]["Waterf_roof"],
    )

    HVAC_Schedule = (;
        Hequip=FT(input_vars["HVACSchedule"]["Hequip"]),
        Hpeople=FT(input_vars["HVACSchedule"]["Hpeople"]),
        LEequip=FT(input_vars["HVACSchedule"]["LEequip"]),
        LEpeople=FT(input_vars["HVACSchedule"]["LEpeople"]),
        AirConRoomFraction=FT(input_vars["HVACSchedule"]["AirConRoomFraction"]),
    )

    TempVec_ittm2Ext = (;
        TWallSun=input_vars["TempVec_ittm2Ext"]["TWallSun"],
        TWallShade=input_vars["TempVec_ittm2Ext"]["TWallShade"],
        TWallIntSun=input_vars["TempVec_ittm2Ext"]["TWallIntSun"],
        TWallIntShade=input_vars["TempVec_ittm2Ext"]["TWallIntShade"],
        TGroundImp=input_vars["TempVec_ittm2Ext"]["TGroundImp"],
        TGroundBare=input_vars["TempVec_ittm2Ext"]["TGroundBare"],
        TGroundVeg=input_vars["TempVec_ittm2Ext"]["TGroundVeg"],
        TTree=input_vars["TempVec_ittm2Ext"]["TTree"],
        TCanyon=input_vars["TempVec_ittm2Ext"]["TCanyon"],
        TRoofVeg=input_vars["TempVec_ittm2Ext"]["TRoofVeg"],
        TRoofIntVeg=input_vars["TempVec_ittm2Ext"]["TRoofIntVeg"],
        TRoofIntImp=input_vars["TempVec_ittm2Ext"]["TRoofIntImp"],
        TRoofImp=input_vars["TempVec_ittm2Ext"]["TRoofImp"],
    )
    Humidity_ittm2Ext = (;
        CanyonSpecific=input_vars["Humidity_ittm2Ext"]["CanyonSpecific"],
    )
    TempVecB_ittm2Ext = (;
        Tbin=input_vars["TempVecB_ittm2Ext"]["Tbin"],
        qbin=input_vars["TempVecB_ittm2Ext"]["qbin"],
        Tinground=input_vars["TempVecB_ittm2Ext"]["Tinground"],
        Tintmass=input_vars["TempVecB_ittm2Ext"]["Tintmass"],
        Tinwallsun=input_vars["TempVecB_ittm2Ext"]["Tinwallsun"],
        Tinwallshd=input_vars["TempVecB_ittm2Ext"]["Tinwallshd"],
        Tceiling=input_vars["TempVecB_ittm2Ext"]["Tceiling"],
        Twindows=input_vars["TempVecB_ittm2Ext"]["Twindows"],
    )
    Meteo_ittm = (;
        SWRin=FT.(input_vars["Meteo_ittm"]["SWRin"]),
        Rain=FT.(input_vars["Meteo_ittm"]["Rain"]),
    )

    T, fval, exitflag = f_solver_tot(
        TempVec_ittm,
        TempVecB_ittm,
        Humidity_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        TempDamp_ittm,
        view_factor,
        Gemeotry_m,
        FractionsGround,
        FractionsRoof,
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
        ParSoilRoof,
        PropOpticalRoof,
        ParThermalRoof,
        ParVegRoof,
        SunPosition,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        PropOpticalIndoors,
        ParHVAC,
        ParThermalBulidFloor,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        TempVec_ittm2Ext,
        Humidity_ittm2Ext,
        TempVecB_ittm2Ext,
        Meteo_ittm,
        Bool(input_vars["RESPreCalc"]),
        FT(input_vars["fconvPreCalc"]),
        FT(input_vars["fconv"]),
        rsRoofPreCalc,
        rsGroundPreCalc,
        rsTreePreCalc,
        HVAC_Schedule,
    )

    @test T ≈ vec(output_vars["T"]) atol=1e-10
    @test fval ≈ vec(output_vars["fval"]) atol=1e-10
    @test exitflag == output_vars["exitflag"]
end
