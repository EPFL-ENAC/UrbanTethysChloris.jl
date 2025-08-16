using Test
using MAT
using UrbanTethysChloris.ConductiveHeat: conductive_heat_flux_green_roof
using ....TestUtils:
    create_height_dependent_vegetation_parameters,
    create_vegetated_soil_parameters,
    create_location_specific_thermal_properties,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "conductive_heat_functions.ConductiveHeatFlux_GreenRoof.mat"
)

# Create parameter structs from input data
ParVegRoof = create_height_dependent_vegetation_parameters(
    FT;
    Rrootl=input_vars["ParVegRoof"]["Rrootl"],
    PsiL50=input_vars["ParVegRoof"]["PsiL50"],
    PsiX50=input_vars["ParVegRoof"]["PsiX50"],
    CASE_ROOT=Int(input_vars["ParVegRoof"]["CASE_ROOT"]),
    ZR95=input_vars["ParVegRoof"]["ZR95"],
    ZR50=input_vars["ParVegRoof"]["ZR50"],
    ZRmax=input_vars["ParVegRoof"]["ZRmax"],
)

ParSoilRoof = create_vegetated_soil_parameters(
    FT;
    Pcla=input_vars["ParSoilRoof"]["Pcla"],
    Psan=input_vars["ParSoilRoof"]["Psan"],
    Porg=input_vars["ParSoilRoof"]["Porg"],
    Kfc=input_vars["ParSoilRoof"]["Kfc"],
    Phy=input_vars["ParSoilRoof"]["Phy"],
    SPAR=Int(input_vars["ParSoilRoof"]["SPAR"]),
    Kbot=input_vars["ParSoilRoof"]["Kbot"],
    Zs=vec(input_vars["ParSoilRoof"]["Zs"]),
    dz1=input_vars["ParSoilRoof"]["dz1"],
    dz2=input_vars["ParSoilRoof"]["dz2"],
)

ParThermalRoof = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalRoof"]["cv_s_imp"],
    lan_dry=input_vars["ParThermalRoof"]["lan_dry_imp"],
)

TempVec_ittm = (; TRoofIntVeg=input_vars["TempVec_ittm"]["TRoofIntVeg"])
Anthropogenic = (; Tb=input_vars["Anthropogenic"]["Tb"])
Owater = (; OwRoofSoilVeg=vec(input_vars["Owater"]["OwRoofSoilVeg"]))
ParCalculation = (; dts=input_vars["ParCalculation"]["dts"])

@testset "MATLAB" begin
    G1, G2, dS = conductive_heat_flux_green_roof(
        vec(input_vars["TemperatureR"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
        Anthropogenic,
        Owater,
        ParVegRoof,
        ParSoilRoof,
        ParThermalRoof,
        ParCalculation,
        Bool(input_vars["BEM_on"]),
    )

    @test G1 ≈ output_vars["G1"]
    @test G2 ≈ output_vars["G2"]
    @test dS ≈ output_vars["dS"]
end
