using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: conductive_heat_flux_building_floor
using ....TestUtils: create_thermal_building, load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "BuildingEnergyModel.ConductiveHeatFluxFR_BuildingFloor.json"
)

ParCalculation = (dts=Int(input_vars["ParCalculation"]["dts"]),)
ParThermalBuildFloor = create_thermal_building(
    FT;
    lan_ground_floor=input_vars["ParThermalBulidFloor"]["lan_ground_floor"],
    cv_ground_floor=input_vars["ParThermalBulidFloor"]["cv_ground_floor"],
)

@testset "MATLAB" begin
    G, Tdp = conductive_heat_flux_building_floor(
        input_vars["Tinground"],
        input_vars["TingroundDamptm1"],
        input_vars["Tingroundtm1"],
        ParCalculation,
        ParThermalBuildFloor,
    )

    @test G ≈ output_vars["G"]
    @test Tdp ≈ output_vars["Tdp"]
end
