using SafeTestsets

@safetestset "water_roof" begin
    include("water_roof.jl")
end

@safetestset "water_vegetation" begin
    include("water_vegetation.jl")
end

@safetestset "water_canyon" begin
    include("water_canyon.jl")
end
