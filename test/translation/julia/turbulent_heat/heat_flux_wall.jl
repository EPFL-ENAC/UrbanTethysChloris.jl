using Test
using MAT
using UrbanTethysChloris.TurbulentHeat: heat_flux_wall
using ....TestUtils:
    create_height_dependent_vegetation_parameters,
    create_location_specific_surface_fractions,
    create_urban_geometry_parameters

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "turbulent_heat_function.HeatFlux_wall.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

Gemeotry_m = create_urban_geometry_parameters(
    FT;
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
    Height_tree=input_vars["Gemeotry_m"]["Height_tree"],
    Radius_tree=input_vars["Gemeotry_m"]["Radius_tree"],
    Hcan_max=input_vars["Gemeotry_m"]["Hcan_max"],
    Hcan_std=input_vars["Gemeotry_m"]["Hcan_std"],
    trees=Bool(input_vars["ParTree"]["trees"]),
)

MeteoData = (;
    Zatm=input_vars["MeteoData"]["Zatm"],
    Tatm=input_vars["MeteoData"]["Tatm"],
    Uatm=input_vars["MeteoData"]["Uatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
)

ParVegTree = create_height_dependent_vegetation_parameters(
    FT; Kopt=input_vars["ParVegTree"]["Kopt"], LAI=input_vars["ParVegTree"]["LAI"]
)

ParVegGround = create_height_dependent_vegetation_parameters(
    FT; hc=input_vars["ParVegGround"]["hc"]
)

FractionsGround = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsGround"]["fveg"],
    fbare=input_vars["FractionsGround"]["fbare"],
    fimp=input_vars["FractionsGround"]["fimp"],
)

@testset "HeatFluxWall" begin
    Hwsun, Hwshade, Ewsun, Ewshade, LEwsun, LEwshade, RES_w1, RES_w2, rap_Zp1_In, rap_Zp2_In, Hwsun1, Hwshade1, Hwsun2, Hwshade2, cp_atm, rho_atm, L_heat, Zp1, Zp2, rap_Zp1, rap_Zp2 = heat_flux_wall(
        vec(input_vars["TemperatureC"]),
        Gemeotry_m,
        MeteoData,
        ParVegTree,
        ParVegGround,
        FractionsGround,
    )

    @test Hwsun ≈ output_vars["Hwsun"]
    @test Hwshade ≈ output_vars["Hwshade"]
    @test Ewsun ≈ output_vars["Ewsun"]
    @test Ewshade ≈ output_vars["Ewshade"]
    @test LEwsun ≈ output_vars["LEwsun"]
    @test LEwshade ≈ output_vars["LEwshade"]
    @test RES_w1 ≈ output_vars["RES_w1"]
    @test RES_w2 ≈ output_vars["RES_w2"]
    @test rap_Zp1_In ≈ output_vars["rap_Zp1_In"]
    @test rap_Zp2_In ≈ output_vars["rap_Zp2_In"]
    @test Hwsun1 ≈ output_vars["Hwsun1"]
    @test Hwshade1 ≈ output_vars["Hwshade1"]
    @test Hwsun2 ≈ output_vars["Hwsun2"]
    @test Hwshade2 ≈ output_vars["Hwshade2"]
    @test cp_atm ≈ output_vars["cp_atm"]
    @test rho_atm ≈ output_vars["rho_atm"]
    @test L_heat ≈ output_vars["L_heat"]
    @test Zp1 ≈ output_vars["Zp1"]
    @test Zp2 ≈ output_vars["Zp2"]
    @test rap_Zp1 ≈ output_vars["rap_Zp1"]
    @test rap_Zp2 ≈ output_vars["rap_Zp2"]
end
