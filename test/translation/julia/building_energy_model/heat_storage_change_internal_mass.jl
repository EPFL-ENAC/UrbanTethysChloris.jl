using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: heat_storage_change_internal_mass
using ....TestUtils:
    create_urban_geometry_parameters, create_thermal_building, load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "BuildingEnergyModel.HeatStorageChangeInternalMass.mat"
)

ParThermalBuildFloor = create_thermal_building(
    FT;
    cv_floor_IntMass=input_vars["ParThermalBuildFloor"]["cv_floor_IntMass"],
    cv_wall_IntMass=input_vars["ParThermalBuildFloor"]["cv_wall_IntMass"],
    dzFloor=input_vars["ParThermalBuildFloor"]["dzFloor"],
    dzWall=input_vars["ParThermalBuildFloor"]["dzWall"],
    FloorHeight=input_vars["ParThermalBuildFloor"]["FloorHeight"],
)

Geometry_m = create_urban_geometry_parameters(
    FT;
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
)

ParCalculation = (dts = Int(input_vars["ParCalculation"]["dts"]))

@testset "MATLAB" begin
    dS = heat_storage_change_internal_mass(
        input_vars["Tintmass"],
        input_vars["Tintmasstm1"],
        ParThermalBuildFloor,
        Geometry_m,
        ParCalculation,
    )

    @test dS ≈ output_vars["dS"]
end
