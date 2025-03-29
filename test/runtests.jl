using SafeTestsets

include("test_utils.jl")

@safetestset "ModelComponents" begin
    include(joinpath("modelcomponents", "ModelComponents.jl"))
end
