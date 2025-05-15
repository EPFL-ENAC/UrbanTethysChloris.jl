using SafeTestsets

@safetestset "evaporation_layers" begin
    include("evaporation_layers.jl")
end

@safetestset "soil_parameters" begin
    include("soil_parameters.jl")
end
