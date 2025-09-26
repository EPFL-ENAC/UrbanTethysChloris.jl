using Test
using MAT
using UrbanTethysChloris.ConductiveHeat: conductive_heat_flux_walls
using ....TestUtils:
    create_location_specific_thermal_properties, create_window_parameters, load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "conductive_heat_functions.ConductiveHeatFlux_Walls.json"
)

# Create parameter structs from input data
ParThermalWall = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalWall"]["cv_s"],
    lan_dry=input_vars["ParThermalWall"]["lan_dry"],
)

ParWindows = create_window_parameters(
    FT;
    Uvalue=input_vars["ParWindows"]["Uvalue"],
    cv_glass=input_vars["ParWindows"]["cv_glass"],
    dztot=input_vars["ParWindows"]["dztot"],
    GlazingRatio=input_vars["ParWindows"]["GlazingRatio"],
)

WallLayers = (;
    dz1_wall=input_vars["WallLayers"]["dz1_wall"],
    dz2_wall=input_vars["WallLayers"]["dz2_wall"],
)

TempVec_ittm = (;
    TWallSun=input_vars["TempVec_ittm"]["TWallSun"],
    TWallShade=input_vars["TempVec_ittm"]["TWallShade"],
    TWallIntSun=input_vars["TempVec_ittm"]["TWallIntSun"],
    TWallIntShade=input_vars["TempVec_ittm"]["TWallIntShade"],
)

TempVecB_ittm = (;
    Tinwallsun=input_vars["TempVecB_ittm"]["Tinwallsun"],
    Tinwallshd=input_vars["TempVecB_ittm"]["Tinwallshd"],
)

Anthropogenic = (; Tb=input_vars["Anthropogenic"]["Tb"])
ParCalculation = (; dts=input_vars["ParCalculation"]["dts"])

@testset "ConductiveHeatFlux_Walls" begin
    G1, G2, dS = conductive_heat_flux_walls(
        vec(input_vars["TemperatureC"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
        TempVecB_ittm,
        Anthropogenic,
        ParThermalWall,
        WallLayers,
        ParCalculation,
        Bool(input_vars["type"]),
        ParWindows,
        Bool(input_vars["BEM_on"]),
    )

    @test G1 ≈ output_vars["G1"]
    @test G2 ≈ output_vars["G2"]
    @test dS ≈ output_vars["dS"]
end
