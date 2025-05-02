using SafeTestsets

@safetestset "AnthropogenicInputs" begin
    include("AnthropogenicInputs.jl")
end

@safetestset "MeteorologicalInputs" begin
    include("MeteorologicalInputs.jl")
end
