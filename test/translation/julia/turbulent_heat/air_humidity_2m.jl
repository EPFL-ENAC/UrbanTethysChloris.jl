using Test
using MAT
using UrbanTethysChloris.TurbulentHeat: air_humidity_2m
using ....TestUtils:
    create_location_specific_surface_fractions,
    create_height_dependent_vegetation_parameters,
    create_urban_geometry_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("turbulent_heat_function.AirHumidity2m.mat")

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "turbulent_heat_function.AirHumidity2m.mat"
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

@testset "AirHumidity2m" begin
    DEi = air_humidity_2m(
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
end
