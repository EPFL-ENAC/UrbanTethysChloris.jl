using Test
using MAT
using UrbanTethysChloris.Water: water_vegetation
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("water_functions.Water_Vegetation.mat")

@testset "Zurich" begin
    q_runon_veg, In_veg, dIn_veg_dt, WBalance_In_veg = water_vegetation(
        input_vars["Rain"],
        input_vars["E_veg"],
        input_vars["In_veg_tm1"],
        input_vars["Sp_In"],
        input_vars["LAI"],
        input_vars["SAI"],
        input_vars["row"],
        input_vars["dth"],
    )

    @test q_runon_veg == output_vars["q_runon_veg"]
    @test In_veg == output_vars["In_veg"]
    @test dIn_veg_dt == output_vars["dIn_veg_dt"]
    @test WBalance_In_veg == output_vars["WBalance_In_veg"]
end
