using SafeTestsets

include("test_utils.jl")

@safetestset "ModelComponents" begin
    include(joinpath("modelcomponents", "ModelComponents.jl"))
end

@safetestset "Radiation" begin
    include(joinpath("radiation", "Radiation.jl"))
end

@safetestset "incoming_longwave" begin
    include("incoming_longwave.jl")
end

@safetestset "set_sun_variables" begin
    include("set_sun_variables.jl")
end
