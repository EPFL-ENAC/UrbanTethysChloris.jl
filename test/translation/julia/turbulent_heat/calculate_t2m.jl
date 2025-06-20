using Test
using MAT
using UrbanTethysChloris.TurbulentHeat: calculate_t2m
using ....TestUtils:
    create_location_specific_surface_fractions,
    create_height_dependent_vegetation_parameters,
    create_urban_geometry_parameters

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "turbulent_heat_function.CalculateT2m.mat"
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

TempVec_ittm = (; T2m=input_vars["TempVec_ittm"]["T2m"])
MeteoData = (; Tatm=input_vars["MeteoData"]["Tatm"])
Gemeotry_m = create_urban_geometry_parameters(
    FT;
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    hcanyon=input_vars["geometry"]["hcanyon"],
)
ParCalculation = (; dts=input_vars["ParCalculation"]["dts"])

@testset "CalculateT2m" begin
    T2m, DHi, Himp_2m, Hbare_2m, Hveg_2m, Hwsun_2m, Hwshade_2m, Hcan_2m = calculate_t2m(
        input_vars["Timp"],
        input_vars["Tbare"],
        input_vars["Tveg"],
        input_vars["Twsun"],
        input_vars["Twshade"],
        input_vars["Tcan"],
        input_vars["Zp1"],
        input_vars["rap_can2m"],
        input_vars["rap_can2m_Inv"],
        input_vars["rb_L"],
        input_vars["RES_w1"],
        FractionsGround,
        Gemeotry_m,
        ParVegGround,
        TempVec_ittm,
        input_vars["cp_atm"],
        input_vars["rho_atm"],
        ParCalculation,
        input_vars["fconv"],
        MeteoData,
    )

    @test T2m ≈ output_vars["T2m"]
    @test DHi ≈ output_vars["DHi"]
    @test Himp_2m ≈ output_vars["Himp_2m"]
    @test Hbare_2m ≈ output_vars["Hbare_2m"]
    @test Hveg_2m ≈ output_vars["Hveg_2m"]
    @test Hwsun_2m ≈ output_vars["Hwsun_2m"]
    @test Hwshade_2m ≈ output_vars["Hwshade_2m"]
    @test Hcan_2m ≈ output_vars["Hcan_2m"]
end
