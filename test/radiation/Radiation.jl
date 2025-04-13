using SafeTestsets

@safetestset "shadow_length" begin
    include("shadow_length.jl")
end

@safetestset "shortwave" begin
    include("shortwave.jl")
end
