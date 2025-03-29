using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    UrbanGeometryParameters, initialize_urbangeometry_parameters

FT = Float64

# Base test case (valid data)
base_data = Dict{String,Any}(
    "Height_canyon" => 6.5,
    "Width_canyon" => 5.78,
    "Width_roof" => 4.73,
    "Radius_tree" => 0.289,
    "Height_tree" => 4.711,
    "Distance_tree" => 0.789,
    "Hcan_max" => NaN,
    "Hcan_std" => NaN,
    "trees" => true,
    "ftree" => 1.0,
)

# Test empty dictionary
@test_throws ArgumentError initialize_urbangeometry_parameters(FT, Dict{String,Any}())

# Test NaN values
for key in [
    "Height_canyon",
    "Width_canyon",
    "Width_roof",
    "Height_tree",
    "Radius_tree",
    "Distance_tree",
]
    nan_data = copy(base_data)
    nan_data[key] = NaN
    @test_throws ArgumentError initialize_urbangeometry_parameters(FT, nan_data)
end

# Test non-positive values
for key in ["Width_canyon", "Width_roof", "Height_canyon"]
    nonsp_data = copy(base_data)
    nonsp_data[key] = 0.0
    @test_throws ArgumentError initialize_urbangeometry_parameters(FT, nonsp_data)
end

# Test invalid tree fraction
invalidftree_data = copy(base_data)
invalidftree_data["ftree"] = 0.5
@test_throws ArgumentError initialize_urbangeometry_parameters(FT, invalidftree_data)

# Test tree-related validations when trees=true
tree_tests = [
    ("Height_tree", 0.0),
    ("Radius_tree", 0.0),
    ("Distance_tree", 0.0),
    ("Radius_tree", 4.0),
    ("Radius_tree", 2.0),
    ("Height_tree", 7.0),
    ("Radius_tree", 5.0),
    ("Radius_tree", 1.0),
    ("Distance_tree", 3.0),
]

for (key, value) in tree_tests
    invalidtree_data = copy(base_data)
    invalidtree_data[key] = value
    @test_throws ArgumentError initialize_urbangeometry_parameters(FT, invalidtree_data)
end

# Test successful initialization
@test_nowarn initialize_urbangeometry_parameters(FT, base_data)

# Test when trees=false
no_trees_data = copy(base_data)
no_trees_data["trees"] = false
@test_nowarn initialize_urbangeometry_parameters(FT, no_trees_data)
