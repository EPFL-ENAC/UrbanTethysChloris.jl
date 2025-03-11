using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    UrbanGeometryParameters, initialize_urbangeometry_parameters

FT = Float64

input_data = Dict{String,Any}(
    "Height_canyon" => 6.5,
    "Width_canyon" => 5.78,
    "Width_roof" => 4.73,
    "Radius_tree" => 100.0,
    "Height_tree" => 0.289,
    "Distance_tree" => 100.0,
    "Hcan_max" => NaN,
    "Hcan_std" => NaN,
    "trees" => true,
    "ftree" => 1.0,
)

ugp = initialize_urbangeometry_parameters(FT, input_data)

@test_throws ArgumentError initialize_urbangeometry_parameters(FT, Dict{String,Any}())
