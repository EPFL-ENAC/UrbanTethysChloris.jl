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
