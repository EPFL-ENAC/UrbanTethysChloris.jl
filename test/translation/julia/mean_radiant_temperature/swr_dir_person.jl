using Test
using MAT
using UrbanTethysChloris.MeanRadiantTemperature: swr_dir_person
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("MRT.SWRDirPerson.mat")

@testset "MATLAB" begin
    SWRdir_Person, SWRdir_in_top, SWRdir_in_bottom, SWRdir_in_east, SWRdir_in_south, SWRdir_in_west, SWRdir_in_north = swr_dir_person(
        input_vars["SWR_dir"],
        input_vars["zeta_S"],
        input_vars["theta_Z"],
        input_vars["BoleanInSun"],
    )

    @test SWRdir_Person ≈ output_vars["SWRdir_Person"]
    @test SWRdir_in_top ≈ output_vars["SWRdir_in_top"]
    @test SWRdir_in_bottom ≈ output_vars["SWRdir_in_bottom"]
    @test SWRdir_in_east ≈ output_vars["SWRdir_in_east"]
    @test SWRdir_in_south ≈ output_vars["SWRdir_in_south"]
    @test SWRdir_in_west ≈ output_vars["SWRdir_in_west"]
    @test SWRdir_in_north ≈ output_vars["SWRdir_in_north"]
end
