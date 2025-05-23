using SafeTestsets

@safetestset "aerodynamic_resistance" begin
    include("aerodynamic_resistance.jl")
end

@safetestset "backcalculate_obhukov_length" begin
    include("backcalculate_obhukov_length.jl")
end

@safetestset "canopy_resistance_an_evolution" begin
    include("canopy_resistance_an_evolution.jl")
end

@safetestset "co2_concentration" begin
    include("co2_concentration.jl")
end

@safetestset "enhancement_factor_ra_pleim" begin
    include("enhancement_factor_ra_pleim.jl")
end

@safetestset "in_canyon_aerodynamic_resistance" begin
    include("in_canyon_aerodynamic_resistance.jl")
end

@safetestset "leaf_boundary_resistance" begin
    include("leaf_boundary_resistance.jl")
end

@safetestset "photosynthesis_biochemical" begin
    include("photosynthesis_biochemical.jl")
end

@safetestset "soil_resistance" begin
    include("soil_resistance.jl")
end
