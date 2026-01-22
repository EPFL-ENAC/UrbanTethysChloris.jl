using Test
using MAT
using UrbanTethysChloris.MeanRadiantTemperature: person_in_shade
using ....TestUtils: create_height_dependent_vegetation_parameters, load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("MRT.PersonInShadeYesOrNo.json")

# Create parameters
ParVegTree = create_height_dependent_vegetation_parameters(
    FT; Kopt=input_vars["ParVegTree"]["Kopt"], LAI=FT(input_vars["ParVegTree"]["LAI"])
)

@testset "MATLAB" begin
    BooleanInSun = person_in_shade(
        Bool(input_vars["trees"]),
        input_vars["h_can"],
        FT(input_vars["w_can"]),
        input_vars["d_tree"],
        input_vars["h_tree"],
        input_vars["r_tree"],
        input_vars["theta_Z"],
        input_vars["theta_n"],
        input_vars["h_P"],
        input_vars["x_P"],
        ParVegTree,
        input_vars["Wcan"],
        input_vars["TimeOfMaxSolAlt"],
        FT(input_vars["TimeHr"]),
    )

    @test BooleanInSun ≈ output_vars["BoleanInSun"]
end
