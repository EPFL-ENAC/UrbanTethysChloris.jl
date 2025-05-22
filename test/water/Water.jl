using SafeTestsets

@safetestset "water_impervious" begin
    include("water_impervious.jl")
end

@safetestset "water_ground" begin
    include("water_ground.jl")
end

@safetestset "water_roof" begin
    include("water_roof.jl")
end

@safetestset "water_soil" begin
    include("water_soil.jl")
end

@safetestset "water_vegetation" begin
    include("water_vegetation.jl")
end
