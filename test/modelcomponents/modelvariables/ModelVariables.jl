using SafeTestsets

@safetestset "BuildingEnergyModelVariables" begin
    include("BuildingEnergyModelVariables.jl")
end

@safetestset "EnvironmentalConditions" begin
    include("EnvironmentalConditions.jl")
end

@safetestset "HeatFluxVariables" begin
    include("HeatFluxVariables.jl")
end

@safetestset "HumidityVariables" begin
    include("HumidityVariables.jl")
end

@safetestset "SolverVariables" begin
    include("SolverVariables.jl")
end

@safetestset "TemperatureVariables" begin
    include("TemperatureVariables.jl")
end
