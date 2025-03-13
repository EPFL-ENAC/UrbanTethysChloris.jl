using SafeTestsets

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

@safetestset "ParameterSet" begin
    include("ParameterSet.jl")
end
