using SafeTestsets

@safetestset "water" begin
    include(joinpath("water", "Water.jl"))
end
