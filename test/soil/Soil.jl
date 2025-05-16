using SafeTestsets

@safetestset "conductivity_suction" begin
    include("conductivity_suction.jl")
end

@safetestset "evaporation_layers" begin
    include("evaporation_layers.jl")
end

@safetestset "root_fraction_general" begin
    include("root_fraction_general.jl")
end

@safetestset "soil_parameters" begin
    include("soil_parameters.jl")
end

@safetestset "soil_parameters2" begin
    include("soil_parameters2.jl")
end

@safetestset "soil_parameters_total" begin
    include("soil_parameters_total.jl")
end
