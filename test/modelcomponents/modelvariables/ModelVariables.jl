using SafeTestsets

@safetestset "BuildingEnergyModelVariables" begin
    include("BuildingEnergyModelVariables.jl")
end

@safetestset "EnergyBalanceVariables" begin
    include("EnergyBalanceVariables.jl")
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

@safetestset "RadiationFluxVariables" begin
    include("RadiationFluxVariables.jl")
end

@safetestset "TemperatureVariables" begin
    include("TemperatureVariables.jl")
end

@safetestset "WaterFluxVariables" begin
    include("WaterFluxVariables.jl")
end

@safetestset "ModelVariableSet" begin
    include("ModelVariableSet.jl")
end
