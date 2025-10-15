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
