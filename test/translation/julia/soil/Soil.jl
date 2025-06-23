using SafeTestsets

@safetestset "Soil_parameters" begin
    include("soil_parameters.jl")
end

@safetestset "SoilParametersTotal" begin
    include("soil_parameters_total.jl")
end
