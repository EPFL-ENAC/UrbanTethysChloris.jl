using SafeTestsets

@safetestset "water_impervious" begin
    include("water_impervious.jl")
end

@safetestset "water_ground" begin
    include("water_ground.jl")
end

@safetestset "water_soil" begin
    include("water_soil.jl")
end
