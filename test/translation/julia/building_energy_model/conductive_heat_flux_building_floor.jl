using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: conductive_heat_flux_building_floor
using ....TestUtils: create_thermal_building

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "BuildingEnergyModel.ConductiveHeatFluxFR_BuildingFloor.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

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
