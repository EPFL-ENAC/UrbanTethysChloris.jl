using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: eb_solver_building_output
using ....TestUtils:
    create_urban_geometry_parameters,
    create_hvac_parameters,
    create_indoor_optical_properties,
    create_thermal_building,
    create_window_parameters

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "BuildingEnergyModel.EBSolver_BuildingOUTPUT.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

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
    RHsetpointCooling=input_vars["ParHVAC"]["RHsetpointCooling"],
    COPAC=input_vars["ParHVAC"]["COPAC"],
    COPHeat=input_vars["ParHVAC"]["COPHeat"],
    ACH=input_vars["ParHVAC"]["ACH"],
    f_ACLatentToQ=input_vars["ParHVAC"]["f_ACLatentToQ"],
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
    FloorHeight=input_vars["ParThermalBulidingInt"]["FloorHeight"],
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
SWRabs_t = (;)
TempDamp_ittm = (; TDampGroundBuild=input_vars["TempDamp_ittm"]["TDampGroundBuild"])

MeteoData = (
    Tatm=input_vars["MeteoData"]["Tatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
)

HVACSchedule = (
    Hequip=input_vars["HVACSchedule"]["Hequip"],
    Hpeople=input_vars["HVACSchedule"]["Hpeople"],
    LEequip=input_vars["HVACSchedule"]["LEequip"],
    LEpeople=input_vars["HVACSchedule"]["LEpeople"],
    AirConRoomFraction=input_vars["HVACSchedule"]["AirConRoomFraction"],
)

@testset "MATLAB" begin
    HbuildInt, LEbuildInt, GbuildInt, SWRabsB, LWRabsB, Tdpfloor, WasteHeat, EnergyUse, HumidityBuilding, ParACHeat, YBuildInt = eb_solver_building_output(
        vec(input_vars["TemperatureC"]),
        vec(input_vars["TemperatureB"]),
        TempVecB_ittm,
        TempVec_ittm,
        Humidity_ittm,
        MeteoData,
        input_vars["SWRinWsun"],
        input_vars["SWRinWshd"],
        input_vars["G2Roof"],
        input_vars["G2WallSun"],
        input_vars["G2WallShade"],
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

    # Test heat fluxes
    @test HbuildInt.HBinRoof ≈ output_vars["HbuildInt"]["HBinRoof"]
    @test HbuildInt.HbinWallSun ≈ output_vars["HbuildInt"]["HbinWallSun"]
    @test HbuildInt.HbinWallshd ≈ output_vars["HbuildInt"]["HbinWallshd"]
    @test HbuildInt.HBinGround ≈ output_vars["HbuildInt"]["HBinGround"]
    @test HbuildInt.HbinIntMass ≈ output_vars["HbuildInt"]["HbinIntMass"]
    @test HbuildInt.HbuildInSurf ≈ output_vars["HbuildInt"]["HbuildInSurf"]
    @test HbuildInt.Hvent ≈ output_vars["HbuildInt"]["Hvent"]
    @test HbuildInt.Hequip ≈ output_vars["HbuildInt"]["Hequip"]
    @test HbuildInt.Hpeople ≈ output_vars["HbuildInt"]["Hpeople"]
    @test HbuildInt.H_AC_Heat ≈ output_vars["HbuildInt"]["H_AC_Heat"]
    @test HbuildInt.dSH_air ≈ output_vars["HbuildInt"]["dSH_air"]

    # Test latent heat fluxes
    @test LEbuildInt.LEvent ≈ output_vars["LEbuildInt"]["LEvent"]
    @test LEbuildInt.LEequip ≈ output_vars["LEbuildInt"]["LEequip"]
    @test LEbuildInt.LEpeople ≈ output_vars["LEbuildInt"]["LEpeople"]
    @test LEbuildInt.LE_AC_Heat ≈ output_vars["LEbuildInt"]["LE_AC_Heat"]
    @test LEbuildInt.dSLE_air ≈ output_vars["LEbuildInt"]["dSLE_air"]

    # Test conductive heat fluxes
    @test GbuildInt.G2Roof ≈ output_vars["GbuildInt"]["G2Roof"]
    @test GbuildInt.G2WallSun ≈ output_vars["GbuildInt"]["G2WallSun"]
    @test GbuildInt.G2WallShade ≈ output_vars["GbuildInt"]["G2WallShade"]
    @test GbuildInt.Gfloor ≈ output_vars["GbuildInt"]["Gfloor"]
    @test GbuildInt.dSinternalMass ≈ output_vars["GbuildInt"]["dSinternalMass"]

    @test SWRabsB.SWRabsCeiling ≈ output_vars["SWRabsB"]["SWRabsCeiling"]
    @test SWRabsB.SWRabsWallsun ≈ output_vars["SWRabsB"]["SWRabsWallsun"]
    @test SWRabsB.SWRabsWallshd ≈ output_vars["SWRabsB"]["SWRabsWallshd"]
    @test SWRabsB.SWRabsGround ≈ output_vars["SWRabsB"]["SWRabsGround"]
    @test SWRabsB.SWRabsInternalMass ≈ output_vars["SWRabsB"]["SWRabsInternalMass"]

    @test LWRabsB.LWRabsCeiling ≈ output_vars["LWRabsB"]["LWRabsCeiling"]
    @test LWRabsB.LWRabsWallsun ≈ output_vars["LWRabsB"]["LWRabsWallsun"]
    @test LWRabsB.LWRabsWallshd ≈ output_vars["LWRabsB"]["LWRabsWallshd"]
    @test LWRabsB.LWRabsGround ≈ output_vars["LWRabsB"]["LWRabsGround"]
    @test LWRabsB.LWRabsInternalMass ≈ output_vars["LWRabsB"]["LWRabsInternalMass"]

    @test Tdpfloor ≈ output_vars["Tdpfloor"]

    @test WasteHeat.SensibleFromAC_Can ≈ output_vars["WasteHeat"]["SensibleFromAC_Can"]
    @test WasteHeat.LatentFromAC_Can ≈ output_vars["WasteHeat"]["LatentFromAC_Can"]
    @test WasteHeat.WaterFromAC_Can ≈ output_vars["WasteHeat"]["WaterFromAC_Can"]
    @test WasteHeat.SensibleFromHeat_Can ≈ output_vars["WasteHeat"]["SensibleFromHeat_Can"]
    @test WasteHeat.LatentFromHeat_Can ≈ output_vars["WasteHeat"]["LatentFromHeat_Can"]
    @test WasteHeat.SensibleFromVent_Can ≈ output_vars["WasteHeat"]["SensibleFromVent_Can"]
    @test WasteHeat.LatentFromVent_Can ≈ output_vars["WasteHeat"]["LatentFromVent_Can"]
    @test WasteHeat.TotAnthInput_URB ≈ output_vars["WasteHeat"]["TotAnthInput_URB"]

    @test EnergyUse.EnergyForAC ≈ output_vars["EnergyUse"]["EnergyForAC"]
    @test EnergyUse.EnergyForAC_H ≈ output_vars["EnergyUse"]["EnergyForAC_H"]
    @test EnergyUse.EnergyForAC_LE ≈ output_vars["EnergyUse"]["EnergyForAC_LE"]
    @test EnergyUse.EnergyForHeating ≈ output_vars["EnergyUse"]["EnergyForHeating"]

    # Test building humidity
    @test HumidityBuilding.qbin ≈ output_vars["HumidityBuilding"]["qbin"]
    @test HumidityBuilding.esatbin ≈ output_vars["HumidityBuilding"]["esatbin"]
    @test HumidityBuilding.ebin ≈ output_vars["HumidityBuilding"]["ebin"]
    @test HumidityBuilding.RHbin ≈ output_vars["HumidityBuilding"]["RHbin"]

    # Test HVAC parameters
    @test ParACHeat.AC_on == Bool(output_vars["ParACHeat"]["AC_on"])
    @test ParACHeat.AC_onCool == Bool(output_vars["ParACHeat"]["AC_onCool"])
    @test ParACHeat.AC_onDehum == Bool(output_vars["ParACHeat"]["AC_onDehum"])
    @test ParACHeat.Heat_on == Bool(output_vars["ParACHeat"]["Heat_on"])

    @test all(isapprox.(YBuildInt, vec(output_vars["YBuildInt"]), atol=1e-11))
end
