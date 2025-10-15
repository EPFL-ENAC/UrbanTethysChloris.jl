using SafeTestsets

@safetestset "EnvironmentalConditions" begin
    include("EnvironmentalConditions.jl")
end

@safetestset "SolverVariables" begin
    include("SolverVariables.jl")
end

@safetestset "TemperatureVariables" begin
    include("TemperatureVariables.jl")
end

@safetestset "Humidity" begin
    include("Humidity.jl")
end

@safetestset "Results2m" begin
    include("Results2m.jl")
end

@safetestset "BuildingEnergyModelVariables" begin
    include("BuildingEnergyModelVariables.jl")
end
