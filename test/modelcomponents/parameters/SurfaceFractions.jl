using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    SurfaceFractions,
    initialize_surfacefractions,
    initialize_locationspecific_surfacefractions

FT = Float64

input_data = Dict{String,Any}()
input_data["roof"] = Dict{String,Any}("fveg" => 0.5, "fimp" => 0.5, "Per_runoff" => 1.0)
input_data["ground"] = Dict{String,Any}(
    "fveg" => 0.545, "fbare" => 0.0, "fimp" => 0.455, "Per_runoff" => 0.9
)

sf = initialize_surfacefractions(FT, input_data)

@test_throws ArgumentError initialize_surfacefractions(FT, Dict{String,Any}())
