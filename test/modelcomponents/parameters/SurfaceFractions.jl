using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    SurfaceFractions, initialize_surfacefractions

FT = Float64

input_data = Dict{String,Any}(
    "fveg_R" => 0.5,
    "fimp_R" => 0.5,
    "Per_runoff_R" => 1.0,
    "fveg_G" => 0.545,
    "fbare_G" => 0.0,
    "fimp_G" => 0.455,
    "Per_runoff_G" => 0.9,
)

sf = initialize_surfacefractions(FT, input_data)

@test_throws ArgumentError initialize_surfacefractions(FT, Dict{String,Any}())
