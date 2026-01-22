using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: ac_heating_turn_on_off
using ....TestUtils:
    create_urban_geometry_parameters, create_hvac_parameters, load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("BuildingEnergyModel.AC_HeatingTurnOnOff.json")

# Create input NamedTuples from MATLAB data
ParHVAC = create_hvac_parameters(
    FT;
    ACon=Bool(input_vars["ParHVAC"]["ACon"]),
    Heatingon=Bool(input_vars["ParHVAC"]["Heatingon"]),
    TsetpointCooling=input_vars["ParHVAC"]["TsetpointCooling"],
    TsetpointHeating=input_vars["ParHVAC"]["TsetpointHeating"],
    RHsetpointCooling=FT(input_vars["ParHVAC"]["RHsetpointCooling"]),
)

TempVecB_ittm = (
    Tbin=input_vars["TempVecB_ittm"]["Tbin"],
    qbin=input_vars["TempVecB_ittm"]["qbin"],
    Tintmass=input_vars["TempVecB_ittm"]["Tintmass"],
    Tinground=input_vars["TempVecB_ittm"]["Tinground"],
    Tceiling=input_vars["TempVecB_ittm"]["Tceiling"],
    Tinwallshd=input_vars["TempVecB_ittm"]["Tinwallshd"],
    Tinwallsun=input_vars["TempVecB_ittm"]["Tinwallsun"],
)

TempVec_ittm = (TCanyon=input_vars["TempVec_ittm"]["TCanyon"],)
Humidity_ittm = (CanyonSpecific=input_vars["Humidity_ittm"]["CanyonSpecific"],)
MeteoData = (Pre=input_vars["MeteoData"]["Pre"],)
Geometry_m = create_urban_geometry_parameters(
    FT;
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
)

@testset "MATLAB" begin
    new_ParHVAC, orig_ParHVAC = ac_heating_turn_on_off(
        ParHVAC,
        TempVecB_ittm,
        TempVec_ittm,
        Humidity_ittm,
        MeteoData,
        Geometry_m,
        Bool(input_vars["BEM_on"]),
    )

    @test new_ParHVAC.ACon == output_vars["ParHVAC"]["ACon"]
    @test new_ParHVAC.AC_onCool == output_vars["ParHVAC"]["AC_onCool"]
    @test new_ParHVAC.AC_onDehum == output_vars["ParHVAC"]["AC_onDehum"]
    @test new_ParHVAC.Heatingon == output_vars["ParHVAC"]["Heatingon"]
    @test new_ParHVAC.TsetpointCooling ≈ output_vars["ParHVAC"]["TsetpointCooling"]
    @test new_ParHVAC.TsetpointHeating ≈ output_vars["ParHVAC"]["TsetpointHeating"]
end
