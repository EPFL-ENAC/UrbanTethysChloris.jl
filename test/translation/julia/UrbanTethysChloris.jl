using SafeTestsets

@safetestset "resistance" begin
    include(joinpath("resistance", "Resistance.jl"))
end

@safetestset "water" begin
    include(joinpath("water", "Water.jl"))
end
