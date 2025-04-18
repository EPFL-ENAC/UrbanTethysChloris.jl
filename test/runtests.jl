using SafeTestsets

include("test_utils.jl")

@safetestset "ModelComponents" begin
    include(joinpath("modelcomponents", "ModelComponents.jl"))
end

@safetestset "Radiation" begin
    include(joinpath("radiation", "Radiation.jl"))
end
