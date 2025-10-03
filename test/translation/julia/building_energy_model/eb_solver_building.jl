using Test
using MAT
using UrbanTethysChloris.Radiation: RadiationFluxes
using UrbanTethysChloris.BuildingEnergyModel: eb_solver_building
using ....TestUtils:
    create_urban_geometry_parameters,
    create_hvac_parameters,
    create_indoor_optical_properties,
    create_thermal_building,
    create_window_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("BuildingEnergyModel.EBSolver_Building.json")

# Create parameter structs from input data
Geometry_m = create_urban_geometry_parameters(
    FT;
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
)

PropOpticalIndoors = create_indoor_optical_properties(
    FT;
    abc=input_vars["PropOpticalIndoors"]["abc"],
    abw=input_vars["PropOpticalIndoors"]["abw"],
    abg=input_vars["PropOpticalIndoors"]["abg"],
    abm=input_vars["PropOpticalIndoors"]["abm"],
    ec=input_vars["PropOpticalIndoors"]["ec"],
    eg=input_vars["PropOpticalIndoors"]["eg"],
    ew=input_vars["PropOpticalIndoors"]["ew"],
    em=input_vars["PropOpticalIndoors"]["em"],
)

ParHVAC = create_hvac_parameters(
    FT;
    ACon=Bool(input_vars["ParHVAC"]["ACon"]),
    Heatingon=Bool(input_vars["ParHVAC"]["Heatingon"]),
    TsetpointCooling=input_vars["ParHVAC"]["TsetpointCooling"],
    TsetpointHeating=input_vars["ParHVAC"]["TsetpointHeating"],
    RHsetpointCooling=FT(input_vars["ParHVAC"]["RHsetpointCooling"]),
    COPAC=input_vars["ParHVAC"]["COPAC"],
    COPHeat=input_vars["ParHVAC"]["COPHeat"],
    ACH=input_vars["ParHVAC"]["ACH"],
    f_ACLatentToQ=FT(input_vars["ParHVAC"]["f_ACLatentToQ"]),
    AC_onCool=Bool(input_vars["ParHVAC"]["AC_onCool"]),
    AC_onDehum=Bool(input_vars["ParHVAC"]["AC_onDehum"]),
)

ParThermalBuildingInt = create_thermal_building(
    FT;
    IntMassOn=Bool(input_vars["ParThermalBulidingInt"]["IntMassOn"]),
    lan_ground_floor=input_vars["ParThermalBulidingInt"]["lan_ground_floor"],
    cv_ground_floor=input_vars["ParThermalBulidingInt"]["cv_ground_floor"],
    cv_floor_IntMass=input_vars["ParThermalBulidingInt"]["cv_floor_IntMass"],
    cv_wall_IntMass=input_vars["ParThermalBulidingInt"]["cv_wall_IntMass"],
    dzFloor=input_vars["ParThermalBulidingInt"]["dzFloor"],
    dzWall=input_vars["ParThermalBulidingInt"]["dzWall"],
    FloorHeight=FT(input_vars["ParThermalBulidingInt"]["FloorHeight"]),
)

ParWindows = create_window_parameters(
    FT; GlazingRatio=input_vars["ParWindows"]["GlazingRatio"]
)

ParCalculation = (
    dts=Int(input_vars["ParCalculation"]["dts"]),
    dth=Int(input_vars["ParCalculation"]["dth"]),
)

TempVecB_ittm = (
    Tbin=input_vars["TempVecB_ittm"]["Tbin"],
    qbin=input_vars["TempVecB_ittm"]["qbin"],
    Tinground=input_vars["TempVecB_ittm"]["Tinground"],
    Tintmass=input_vars["TempVecB_ittm"]["Tintmass"],
    Tinwallsun=input_vars["TempVecB_ittm"]["Tinwallsun"],
    Tinwallshd=input_vars["TempVecB_ittm"]["Tinwallshd"],
    Tceiling=input_vars["TempVecB_ittm"]["Tceiling"],
)

TempVec_ittm = (TWallIntSun=input_vars["TempVec_ittm"]["TWallIntSun"],)
Humidity_ittm = (;)
SWRabs_t = RadiationFluxes{FT}(;
    GroundImp=input_vars["SWRabs_t"]["SWRabsGroundImp"],
    GroundBare=input_vars["SWRabs_t"]["SWRabsGroundBare"],
    GroundVeg=input_vars["SWRabs_t"]["SWRabsGroundVeg"],
    Tree=input_vars["SWRabs_t"]["SWRabsTree"],
    WallSun=input_vars["SWRabs_t"]["SWRabsWallSun"],
    WallShade=input_vars["SWRabs_t"]["SWRabsWallShade"],
    TotalGround=input_vars["SWRabs_t"]["SWRabsTotalGround"],
    TotalCanyon=input_vars["SWRabs_t"]["SWRabsTotalCanyon"],
)
TempDamp_ittm = (; TDampGroundBuild=input_vars["TempDamp_ittm"]["TDampGroundBuild"])

MeteoData = (
    Tatm=input_vars["MeteoData"]["Tatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
)

HVACSchedule = (
    Hequip=FT(input_vars["HVACSchedule"]["Hequip"]),
    Hpeople=FT(input_vars["HVACSchedule"]["Hpeople"]),
    LEequip=FT(input_vars["HVACSchedule"]["LEequip"]),
    LEpeople=FT(input_vars["HVACSchedule"]["LEpeople"]),
    AirConRoomFraction=input_vars["HVACSchedule"]["AirConRoomFraction"],
)

@testset "MATLAB" begin
    YBuildInt, WasteHeat = eb_solver_building(
        vec(input_vars["TemperatureC"]),
        vec(input_vars["TemperatureB"]),
        TempVecB_ittm,
        TempVec_ittm,
        Humidity_ittm,
        MeteoData,
        FT(input_vars["SWRinWsun"]),
        FT(input_vars["SWRinWshd"]),
        FT(input_vars["G2Roof"]),
        FT(input_vars["G2WallSun"]),
        FT(input_vars["G2WallShade"]),
        TempDamp_ittm,
        SWRabs_t,
        Geometry_m,
        PropOpticalIndoors,
        ParHVAC,
        ParCalculation,
        ParThermalBuildingInt,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        HVACSchedule,
    )

    @test all(isapprox.(YBuildInt, vec(output_vars["YBuildInt"]), atol=1e-11))
    @test WasteHeat.SensibleFromAC_Can ≈ output_vars["WasteHeat"]["SensibleFromAC_Can"]
    @test WasteHeat.LatentFromAC_Can ≈ output_vars["WasteHeat"]["LatentFromAC_Can"]
    @test WasteHeat.WaterFromAC_Can ≈ output_vars["WasteHeat"]["WaterFromAC_Can"]
    @test WasteHeat.SensibleFromHeat_Can ≈ output_vars["WasteHeat"]["SensibleFromHeat_Can"]
    @test WasteHeat.LatentFromHeat_Can ≈ output_vars["WasteHeat"]["LatentFromHeat_Can"]
    @test WasteHeat.SensibleFromVent_Can ≈ output_vars["WasteHeat"]["SensibleFromVent_Can"]
    @test WasteHeat.LatentFromVent_Can ≈ output_vars["WasteHeat"]["LatentFromVent_Can"]
    @test WasteHeat.TotAnthInput_URB ≈ output_vars["WasteHeat"]["TotAnthInput_URB"]
end
