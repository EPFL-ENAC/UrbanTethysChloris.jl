using Test
using UrbanTethysChloris.Soil: soil_moistures_rich_comp_lat3
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "soil_functions.SOIL_MOISTURES_RICH_COMP_LAT3.json"
)

FT = Float64

@testset "MATLAB" begin
    alpVG = FT(input_vars["alpVG"])
    Cimp = FT(input_vars["Cimp"])
    Cbare = FT(input_vars["Cbare"])
    Cveg = FT(input_vars["Cveg"])
    dz = FT(input_vars["dz"])
    fbare = FT(input_vars["fbare"])
    fimp = FT(input_vars["fimp"])
    fveg = FT(input_vars["fveg"])
    Ks = FT(input_vars["Ks"])
    L = FT(input_vars["L"])
    nVG = FT(input_vars["nVG"])
    O33 = input_vars["O33"]
    Ohy = input_vars["Ohy"]
    Osat = input_vars["Osat"]
    Pe = input_vars["Pe"]
    Vlat = input_vars["Vlat"]
    Wcan = input_vars["Wcan"]

    @testset "SPAR 2" begin
        dVlat = soil_moistures_rich_comp_lat3(
            Vlat,
            dz,
            2,
            Ks,
            Osat,
            Ohy,
            L,
            Pe,
            O33,
            alpVG,
            nVG,
            Cimp,
            Cbare,
            Cveg,
            fimp,
            fbare,
            fveg,
            Wcan,
        )

        @test dVlat ≈ output_vars["dVlat"]
    end
end
