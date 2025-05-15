using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    LocationProperties, initialize_locationproperties

FT = Float64

location_data = Dict{String,Any}(
    "phi" => 25.3, "lambda" => 55.4, "theta_canyon" => deg2rad(45.0), "DeltaGMT" => 4.0
)

@test_nowarn initialize_locationproperties(FT, location_data)
