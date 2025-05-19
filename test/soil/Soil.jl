using SafeTestsets

@safetestset "conductivity_suction" begin
    include("conductivity_suction.jl")
end

@safetestset "evaporation_layers" begin
    include("evaporation_layers.jl")
end

@safetestset "infiltration" begin
    include("infiltration.jl")
end

@safetestset "leakage_bottom" begin
    include("leakage_bottom.jl")
end

@safetestset "root_fraction_general" begin
    include("root_fraction_general.jl")
end

@safetestset "root_soil_conductance" begin
    include("root_soil_conductance.jl")
end

@safetestset "soil_moistures_rich_comp" begin
    include("soil_moistures_rich_comp.jl")
end

@safetestset "soil_moistures_rich_comp_lat2" begin
    include("soil_moistures_rich_comp_lat2.jl")
end

@safetestset "soil_moistures_rich_comp_lat3" begin
    include("soil_moistures_rich_comp_lat3.jl")
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

@safetestset "soil_thermal_properties" begin
    include("soil_thermal_properties.jl")
end

@safetestset "soil_water_multilayer" begin
    include("soil_water_multilayer.jl")
end

@safetestset "volume_correction" begin
    include("volume_correction.jl")
end
