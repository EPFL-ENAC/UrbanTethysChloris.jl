using Test
using MAT
using UrbanTethysChloris.ConductiveHeat: conductive_heat_flux_ground_vb
using ....TestUtils:
    create_location_specific_thermal_properties,
    create_location_specific_surface_fractions,
    create_vegetated_soil_parameters,
    create_height_dependent_vegetation_parameters

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "conductive_heat_functions.ConductiveHeatFluxFR_GroundVegBare.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

# Create parameter structs from input data
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
    ZR95=[input_vars["ParVegTree"]["ZR95"]],
    ZR50=[input_vars["ParVegTree"]["ZR50"]],
    ZRmax=[input_vars["ParVegTree"]["ZRmax"]],
)

ParVegGround = create_height_dependent_vegetation_parameters(
    FT;
    CASE_ROOT=Int(input_vars["ParVegGround"]["CASE_ROOT"]),
    ZR95=[input_vars["ParVegGround"]["ZR95"]],
    ZR50=[input_vars["ParVegGround"]["ZR50"]],
    ZRmax=[input_vars["ParVegGround"]["ZRmax"]],
)

TempDamp_ittm = (;
    TDampGroundBare=input_vars["TempDamp_ittm"]["TDampGroundBare"],
    TDampGroundVeg=input_vars["TempDamp_ittm"]["TDampGroundVeg"],
)
TempVec_ittm = (;
    TGroundBare=input_vars["TempVec_ittm"]["TGroundBare"],
    TGroundVeg=input_vars["TempVec_ittm"]["TGroundVeg"],
)
Owater_ittm = (;
    OwGroundSoilBare=vec(input_vars["Owater_ittm"]["OwGroundSoilBare"]),
    OwGroundSoilVeg=vec(input_vars["Owater_ittm"]["OwGroundSoilVeg"]),
)
ParCalculation = (; dts=Int(input_vars["ParCalculation"]["dts"]))

@testset "ConductiveHeatFluxFR_GroundVegBare" begin
    G, Tdp = conductive_heat_flux_ground_vb(
        vec(input_vars["TemperatureC"]),
        TempDamp_ittm,
        Owater_ittm,
        TempVec_ittm,
        ParCalculation,
        ParSoilGround,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        Int(input_vars["type"]),
    )

    @test G ≈ output_vars["G"]
    @test Tdp ≈ output_vars["Tdp"]
end
