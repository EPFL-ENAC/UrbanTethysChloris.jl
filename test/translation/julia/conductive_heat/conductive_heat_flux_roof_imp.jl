using Test
using MAT
using UrbanTethysChloris.ConductiveHeat: conductive_heat_flux_roof_imp
using ....TestUtils:
    create_vegetated_soil_parameters,
    create_location_specific_thermal_properties,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "conductive_heat_functions.ConductiveHeatFlux_RoofImp.json"
)

# Create parameter structs from input data
ParThermalRoof = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalRoof"]["cv_s_imp"],
    lan_dry=input_vars["ParThermalRoof"]["lan_dry_imp"],
)

ParSoilRoof = create_vegetated_soil_parameters(
    FT; dz1=input_vars["ParSoilRoof"]["dz1"], dz2=input_vars["ParSoilRoof"]["dz2"]
)

TempVec_ittm = (; TRoofIntImp=input_vars["TempVec_ittm"]["TRoofIntImp"])
Anthropogenic = (; Tb=input_vars["Anthropogenic"]["Tb"])
ParCalculation = (; dts=input_vars["ParCalculation"]["dts"])

@testset "MATLAB" begin
    G1, G2, dS = conductive_heat_flux_roof_imp(
        vec(input_vars["TemperatureR"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
        Anthropogenic,
        ParThermalRoof,
        ParSoilRoof,
        ParCalculation,
        Bool(input_vars["BEM_on"]),
    )

    @test G1 ≈ output_vars["G1"]
    @test G2 ≈ output_vars["G2"]
    @test dS ≈ output_vars["dS"]
end
