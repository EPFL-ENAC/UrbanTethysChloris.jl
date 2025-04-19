using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    initialize_urbangeometry_parameters, initialize_person_parameters
using UrbanTethysChloris.RayTracing: view_factors_canyon

FT = Float64

@testset "view_factors_canyon" begin
    @testset "without trees" begin
        geom_data = Dict{String,Any}(
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

        geom = initialize_urbangeometry_parameters(FT, geom_data)

        person_data = Dict{String,Any}(
            "PositionPx" => 2.9,
            "PositionPz" => 1.1,
            "PersonWidth" => 0.14,
            "PersonHeight" => 0.87,
            "HeightWind" => 1.1,
        )
        person = initialize_person_parameters(FT, person_data)

        vf, vfp = view_factors_canyon(geom, person)
    end
end
