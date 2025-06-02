using SafeTestsets

include("test_utils.jl")

@safetestset "ModelComponents" begin
    include(joinpath("modelcomponents", "ModelComponents.jl"))
end

@safetestset "Radiation" begin
    include(joinpath("radiation", "Radiation.jl"))
end

@safetestset "RayTracing" begin
    include(joinpath("ray_tracing", "RayTracing.jl"))
end

@safetestset "Soil" begin
    include(joinpath("soil", "Soil.jl"))
end

@safetestset "Water" begin
    include(joinpath("water", "Water.jl"))
end

@safetestset "incoming_longwave" begin
    include("incoming_longwave.jl")
end

@safetestset "set_sun_variables" begin
    include("set_sun_variables.jl")
end

@safetestset "Translation" begin
    include(joinpath("translation", "julia", "UrbanTethysChloris.jl"))
end
