using Test
using MAT
using UrbanTethysChloris.ConductiveHeat: conductive_heat_flux_ground_fr
using ....TestUtils:
    create_location_specific_thermal_properties,
    create_location_specific_surface_fractions,
    create_vegetated_soil_parameters,
    create_height_dependent_vegetation_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "conductive_heat_functions.ConductiveHeatFluxFR_GroundImp.mat"
)

# Create parameter structs from input data
ParThermalGround = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalGround"]["cv_s_imp"],
    lan_dry=input_vars["ParThermalGround"]["lan_dry_imp"],
)

FractionsGround = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsGround"]["fveg"],
    fbare=input_vars["FractionsGround"]["fbare"],
    fimp=input_vars["FractionsGround"]["fimp"],
)

ParSoilGround = create_vegetated_soil_parameters(
    FT;
    Pcla=input_vars["ParSoilGround"]["Pcla"],
    Psan=input_vars["ParSoilGround"]["Psan"],
    Porg=input_vars["ParSoilGround"]["Porg"],
    Kfc=input_vars["ParSoilGround"]["Kfc"],
    Phy=input_vars["ParSoilGround"]["Phy"],
    SPAR=Int(input_vars["ParSoilGround"]["SPAR"]),
    Kbot=input_vars["ParSoilGround"]["Kbot"],
    Zs=vec(input_vars["ParSoilGround"]["Zs"]),
)

ParVegTree = create_height_dependent_vegetation_parameters(
    FT;
    CASE_ROOT=Int(input_vars["ParVegTree"]["CASE_ROOT"]),
    ZR95=input_vars["ParVegTree"]["ZR95"],
    ZR50=input_vars["ParVegTree"]["ZR50"],
    ZRmax=input_vars["ParVegTree"]["ZRmax"],
)

ParVegGround = create_height_dependent_vegetation_parameters(
    FT;
    CASE_ROOT=Int(input_vars["ParVegGround"]["CASE_ROOT"]),
    ZR95=input_vars["ParVegGround"]["ZR95"],
    ZR50=input_vars["ParVegGround"]["ZR50"],
    ZRmax=input_vars["ParVegGround"]["ZRmax"],
)

TempDamp_ittm = (; TDampGroundImp=input_vars["TempDamp_ittm"]["TDampGroundImp"])
TempVec_ittm = (; TGroundImp=input_vars["TempVec_ittm"]["TGroundImp"])
Owater_ittm = (; OwGroundSoilImp=vec(input_vars["Owater_ittm"]["OwGroundSoilImp"]))
ParCalculation = (; dts=Int(input_vars["ParCalculation"]["dts"]))

@testset "ConductiveHeatFluxFR_GroundImp" begin
    G, Tdp = conductive_heat_flux_ground_fr(
        vec(input_vars["TemperatureC"]),
        TempDamp_ittm,
        TempVec_ittm,
        Owater_ittm,
        ParCalculation,
        ParThermalGround,
        FractionsGround,
        ParSoilGround,
        ParVegTree,
        ParVegGround,
    )

    @test G ≈ output_vars["G"]
    @test Tdp ≈ output_vars["Tdp"]
end
