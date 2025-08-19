using Test
using MAT
using UrbanTethysChloris.TurbulentHeat: heat_flux_canyon
using ....TestUtils:
    create_height_dependent_vegetation_parameters,
    create_urban_geometry_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("turbulent_heat_function.HeatFlux_canyon.mat")

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
    q_atm=input_vars["MeteoData"]["q_atm"],
    ea=input_vars["MeteoData"]["ea"],
)

ParVegTree = create_height_dependent_vegetation_parameters(
    FT; Kopt=input_vars["ParVegTree"]["Kopt"], LAI=input_vars["ParVegTree"]["LAI"]
)

@testset "HeatFluxCanyon" begin
    Hcanyon, LEcanyon, ra_canyon, ra_orig, fconv, HumidityCan = heat_flux_canyon(
        vec(input_vars["TemperatureC"]),
        Gemeotry_m,
        MeteoData,
        ParVegTree,
        input_vars["fconvPreCalc"],
        input_vars["fconv"],
    )

    @test Hcanyon ≈ output_vars["Hcanyon"]
    @test LEcanyon ≈ output_vars["LEcanyon"]
    @test ra_canyon ≈ output_vars["ra_canyon"]
    @test ra_orig ≈ output_vars["ra_orig"]
    @test fconv ≈ output_vars["fconv"]
    @test HumidityCan.CanyonRelative ≈ output_vars["HumidityCan"]["CanyonRelative"]
    @test HumidityCan.CanyonSpecific ≈ output_vars["HumidityCan"]["CanyonSpecific"]
    @test HumidityCan.CanyonVapourPre ≈ output_vars["HumidityCan"]["CanyonVapourPre"]
    @test HumidityCan.CanyonRelativeSat ≈ output_vars["HumidityCan"]["CanyonRelativeSat"]
    @test HumidityCan.CanyonSpecificSat ≈ output_vars["HumidityCan"]["CanyonSpecificSat"]
    @test HumidityCan.CanyonVapourPreSat ≈ output_vars["HumidityCan"]["CanyonVapourPreSat"]
end
