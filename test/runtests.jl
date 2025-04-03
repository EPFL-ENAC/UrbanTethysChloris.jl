using SafeTestsets

include("test_utils.jl")

@safetestset "utils" begin
    include("utils.jl")
end

@safetestset "ModelComponents" begin
    include(joinpath("modelcomponents", "ModelComponents.jl"))
end
