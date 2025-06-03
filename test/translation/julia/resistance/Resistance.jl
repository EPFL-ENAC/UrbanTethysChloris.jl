using SafeTestsets

@safetestset "precalculate_stomatal_resistance_roof" begin
    include("precalculate_stomatal_resistance_roof.jl")
end

@safetestset "precalculate_stomatal_resistance_ground_tree" begin
    include("precalculate_stomatal_resistance_ground_tree.jl")
end
