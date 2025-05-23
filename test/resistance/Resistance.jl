using SafeTestsets

@safetestset "aerodynamic_resistance" begin
    include("aerodynamic_resistance.jl")
end

@safetestset "leaf_boundary_resistance" begin
    include("leaf_boundary_resistance.jl")
end

@safetestset "photosynthesis_biochemical" begin
    include("photosynthesis_biochemical.jl")
end
