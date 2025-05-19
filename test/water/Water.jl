using SafeTestsets

@safetestset "water_impervious" begin
    include("water_impervious.jl")
end

@safetestset "water_ground" begin
    include("water_ground.jl")
end
