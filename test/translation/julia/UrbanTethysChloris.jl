using SafeTestsets

@safetestset "BuildingEnergyModel" begin
    include(joinpath("building_energy_model", "BuildingEnergyModel.jl"))
end

@safetestset "conductive_heat" begin
    include(joinpath("conductive_heat", "ConductiveHeat.jl"))
end

@safetestset "mean_radiant_temperature" begin
    include(joinpath("mean_radiant_temperature", "MeanRadiantTemperature.jl"))
end

@safetestset "outdoor_thermal_comfort" begin
    include(joinpath("outdoor_thermal_comfort", "OutdoorThermalComfort.jl"))
end

@safetestset "ray_tracing" begin
    include(joinpath("ray_tracing", "RayTracing.jl"))
end

@safetestset "resistance" begin
    include(joinpath("resistance", "Resistance.jl"))
end

@safetestset "soil" begin
    include(joinpath("soil", "Soil.jl"))
end

@safetestset "turbulent_heat" begin
    include(joinpath("turbulent_heat", "TurbulentHeat.jl"))
end

@safetestset "water" begin
    include(joinpath("water", "Water.jl"))
end

@safetestset "eb_solver_roof" begin
    include("eb_solver_roof.jl")
end

@safetestset "eb_solver_canyon" begin
    include("eb_solver_canyon.jl")
end

@safetestset "eb_wb_roof" begin
    include("eb_wb_roof.jl")
end

@safetestset "eb_solver_urban_climate_building_energy_model" begin
    include("eb_solver_urban_climate_building_energy_model.jl")
end

@safetestset "f_solver_tot" begin
    include("f_solver_tot.jl")
end

# TODO: remove ./julia subdirectory level, go straigth for test/transation/UrbanTethysChloris.jl
