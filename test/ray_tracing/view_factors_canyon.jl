using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    initialize_urbangeometry_parameters, initialize_person_parameters
using UrbanTethysChloris.RayTracing: view_factors_canyon
using Random

FT = Float64

@testset "view_factors_canyon" begin
    person_data = Dict{String,Any}(
        "PositionPx" => 2.9,
        "PositionPz" => 1.1,
        "PersonWidth" => 0.14,
        "PersonHeight" => 0.87,
        "HeightWind" => 1.1,
    )
    person = initialize_person_parameters(FT, person_data)

    geom_data = Dict{String,Any}(
        "Height_canyon" => 6.5,
        "Width_canyon" => 5.78,
        "Width_roof" => 4.73,
        "Radius_tree" => 0.289,
        "Height_tree" => 4.711,
        "Distance_tree" => 0.789,
        "Hcan_max" => NaN,
        "Hcan_std" => NaN,
        "ftree" => 1.0,
    )
    @testset "with trees" begin
        geom_data["trees"] = true
        geom = initialize_urbangeometry_parameters(FT, geom_data)

        vf, vfp = view_factors_canyon(geom, person)
    end

    @testset "without trees" begin
        geom_data["trees"] = false
        geom = initialize_urbangeometry_parameters(FT, geom_data)

        Random.seed!(123)  # Set seed for reproducibility
        vf, vfp = view_factors_canyon(geom, person)
    end
end
