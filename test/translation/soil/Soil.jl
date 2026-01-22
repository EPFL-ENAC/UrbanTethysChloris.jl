using SafeTestsets

@safetestset "SOIL_MOISTURES_RICH_COMP_LAT3" begin
    include("soil_moistures_rich_comp_lat3.jl")
end

@safetestset "Soil_parameters" begin
    include("soil_parameters.jl")
end

@safetestset "SoilParametersTotal" begin
    include("soil_parameters_total.jl")
end
