using SafeTestsets

@safetestset "SoilParameters.jl" begin
    include("SoilParameters.jl")
end

@safetestset "SurfaceFractions" begin
    include("SurfaceFractions.jl")
end

@safetestset "UrbanGeometryParameters" begin
    include("UrbanGeometryParameters.jl")
end

@safetestset "ThermalProperties" begin
    include("ThermalProperties.jl")
end

@safetestset "VegetationParameters" begin
    include("VegetationParameters.jl")
end

@safetestset "OpticalProperties" begin
    include("OpticalProperties.jl")
end

@safetestset "PersonParameters" begin
    include("PersonParameters.jl")
end

@safetestset "LocationProperties" begin
    include("LocationProperties.jl")
end

@testset "Parameters" begin
    include("BuildingEnergyModelParameters.jl")
end

@safetestset "ParameterSet" begin
    include("ParameterSet.jl")
end
