using Test
using MAT
using UrbanTethysChloris.TurbulentHeat: air_humidity_2m_output
using ....TestUtils:
    create_location_specific_surface_fractions,
    create_height_dependent_vegetation_parameters,
    create_urban_geometry_parameters

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "turbulent_heat_function.AirHumidity2mOutput.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

FractionsGround = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsGround"]["fveg"],
    fbare=input_vars["FractionsGround"]["fbare"],
    fimp=input_vars["FractionsGround"]["fimp"],
)

ParVegGround = create_height_dependent_vegetation_parameters(
    FT; LAI=input_vars["ParVegGround"]["LAI"], SAI=input_vars["ParVegGround"]["SAI"]
)

Humidity_ittm = (; q2m=input_vars["Humidity_ittm"]["q2m"])
MeteoData = (; Tatm=input_vars["MeteoData"]["Tatm"])
Gemeotry_m = create_urban_geometry_parameters(
    FT;
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
)
ParCalculation = (; dts=input_vars["ParCalculation"]["dts"])

@testset "AirHumidity2mOutput" begin
    DEi, Eimp_2m, Ebare_soil_2m, Eveg_int_2m, Eveg_soil_2m, TEveg_2m, Ecan_2m, q2m, e_T2m, RH_T2m, qcan, e_Tcan, RH_Tcan = air_humidity_2m_output(
        input_vars["q2m"],
        input_vars["T2m"],
        input_vars["Timp"],
        input_vars["Tbare"],
        input_vars["Tveg"],
        input_vars["Tcan"],
        input_vars["qcan"],
        input_vars["rap_can2m"],
        input_vars["rap_can2m_Inv"],
        input_vars["rb_L"],
        input_vars["alp_soil_bare"],
        input_vars["r_soil_bare"],
        input_vars["alp_soil_veg"],
        input_vars["r_soil_veg"],
        input_vars["rs_sun_L"],
        input_vars["rs_shd_L"],
        input_vars["dw_L"],
        input_vars["Fsun_L"],
        input_vars["Fshd_L"],
        FractionsGround,
        ParVegGround,
        input_vars["Eimp"],
        input_vars["Ebare"],
        input_vars["Eveg_int"],
        input_vars["Eveg_pond"],
        input_vars["Eveg_soil"],
        input_vars["TEveg"],
        input_vars["Pre"],
        Humidity_ittm,
        input_vars["fconv"],
        MeteoData,
        Gemeotry_m,
        input_vars["rho_atm"],
        input_vars["Zp1"],
        ParCalculation,
    )

    @test DEi ≈ output_vars["DEi"]
    @test Eimp_2m ≈ output_vars["Eimp_2m"]
    @test Ebare_soil_2m ≈ output_vars["Ebare_soil_2m"]
    @test Eveg_int_2m ≈ output_vars["Eveg_int_2m"]
    @test Eveg_soil_2m ≈ output_vars["Eveg_soil_2m"]
    @test TEveg_2m ≈ output_vars["TEveg_2m"]
    @test Ecan_2m ≈ output_vars["Ecan_2m"]
    @test q2m ≈ output_vars["q2m"]
    @test e_T2m ≈ output_vars["e_T2m"]
    @test RH_T2m ≈ output_vars["RH_T2m"]
    @test qcan ≈ output_vars["qcan"]
    @test e_Tcan ≈ output_vars["e_Tcan"]
    @test RH_Tcan ≈ output_vars["RH_Tcan"]
end
