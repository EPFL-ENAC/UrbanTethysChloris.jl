using Test
using MAT
using UrbanTethysChloris.TurbulentHeat: heat_flux_roof
using ....TestUtils:
    create_height_dependent_vegetation_parameters,
    create_location_specific_surface_fractions,
    create_urban_geometry_parameters,
    create_vegetated_soil_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("turbulent_heat_function.HeatFlux_roof.mat")

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "turbulent_heat_function.HeatFlux_roof.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

# Create structured inputs with correct types
FractionsRoof = create_location_specific_surface_fractions(
    FT; fveg=input_vars["FractionsRoof"]["fveg"], fimp=input_vars["FractionsRoof"]["fimp"]
)

ParVegRoof = create_height_dependent_vegetation_parameters(
    FT;
    LAI=input_vars["ParVegRoof"]["LAI"],
    SAI=input_vars["ParVegRoof"]["SAI"],
    hc=input_vars["ParVegRoof"]["hc"],
    d_leaf=input_vars["ParVegRoof"]["d_leaf"],
    Kopt=input_vars["ParVegRoof"]["Kopt"],
    Knit=input_vars["ParVegRoof"]["Knit"],
    Psi_sto_50=input_vars["ParVegRoof"]["Psi_sto_50"],
    Psi_sto_00=input_vars["ParVegRoof"]["Psi_sto_00"],
    CT=Int(input_vars["ParVegRoof"]["CT"]),
    Vmax=input_vars["ParVegRoof"]["Vmax"],
    DSE=input_vars["ParVegRoof"]["DSE"],
    Ha=input_vars["ParVegRoof"]["Ha"],
    FI=input_vars["ParVegRoof"]["FI"],
    Do=input_vars["ParVegRoof"]["Do"],
    a1=input_vars["ParVegRoof"]["a1"],
    go=input_vars["ParVegRoof"]["go"],
    e_rel=input_vars["ParVegRoof"]["e_rel"],
    e_relN=input_vars["ParVegRoof"]["e_relN"],
    gmes=input_vars["ParVegRoof"]["gmes"],
    rjv=input_vars["ParVegRoof"]["rjv"],
    mSl=input_vars["ParVegRoof"]["mSl"],
    Sl=input_vars["ParVegRoof"]["Sl"],
    CASE_ROOT=Int(input_vars["ParVegRoof"]["CASE_ROOT"]),
    ZR95=input_vars["ParVegRoof"]["ZR95"],
    ZR50=input_vars["ParVegRoof"]["ZR50"],
    ZRmax=input_vars["ParVegRoof"]["ZRmax"],
)

ParSoilRoof = create_vegetated_soil_parameters(
    FT;
    Pcla=input_vars["ParSoilRoof"]["Pcla"],
    Psan=input_vars["ParSoilRoof"]["Psan"],
    Porg=input_vars["ParSoilRoof"]["Porg"],
    Kfc=input_vars["ParSoilRoof"]["Kfc"],
    Phy=input_vars["ParSoilRoof"]["Phy"],
    SPAR=Int(input_vars["ParSoilRoof"]["SPAR"]),
    Kbot=input_vars["ParSoilRoof"]["Kbot"],
    Sp_In=input_vars["ParSoilRoof"]["Sp_In"],
    Zs=vec(input_vars["ParSoilRoof"]["Zs"]),
)

Gemeotry_m = create_urban_geometry_parameters(
    FT; Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"]
)

MeteoData = (;
    Zatm=input_vars["MeteoData"]["Zatm"],
    Tatm=input_vars["MeteoData"]["Tatm"],
    Uatm=input_vars["MeteoData"]["Uatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
    Catm_O2=input_vars["MeteoData"]["Catm_O2"],
    Catm_CO2=input_vars["MeteoData"]["Catm_CO2"],
)

TempVec_ittm = (; TRoofVeg=input_vars["TempVec_ittm"]["TRoofVeg"])
HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])
ParCalculation = (;
    dth=input_vars["ParCalculation"]["dth"], row=input_vars["ParCalculation"]["row"]
)

SoilPotW_ittm = (; SoilPotWRoofVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWRoofVeg_L"])
Owater_ittm = (; OwRoofSoilVeg=vec(input_vars["Owater_ittm"]["OwRoofSoilVeg"]))
Vwater_ittm = (; VRoofSoilVeg=vec(input_vars["Vwater_ittm"]["VRoofSoilVeg"]))
ExWater_ittm = (; ExWaterRoofVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterRoofVeg_L"]))
Int_ittm = (;
    IntRoofImp=input_vars["Int_ittm"]["IntRoofImp"],
    IntRoofVegPlant=input_vars["Int_ittm"]["IntRoofVegPlant"],
    IntRoofVegGround=input_vars["Int_ittm"]["IntRoofVegGround"],
)
CiCO2Leaf_ittm = (;
    CiCO2LeafRoofVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegSun"],
    CiCO2LeafRoofVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegShd"],
)

rsRoofPreCalc = (;
    rs_sun=input_vars["rsRoofPreCalc"]["rs_sun"],
    rs_shd=input_vars["rsRoofPreCalc"]["rs_shd"],
    Ci_sun=input_vars["rsRoofPreCalc"]["Ci_sun"],
    Ci_shd=input_vars["rsRoofPreCalc"]["Ci_shd"],
)

@testset "HeatFluxRoof" begin
    results = heat_flux_roof(
        vec(input_vars["TemperatureR"]),
        TempVec_ittm,
        MeteoData,
        HumidityAtm,
        ParVegRoof,
        FractionsRoof,
        Gemeotry_m,
        ParSoilRoof,
        ParCalculation,
        SoilPotW_ittm,
        Owater_ittm,
        Vwater_ittm,
        ExWater_ittm,
        Int_ittm,
        CiCO2Leaf_ittm,
        FT(input_vars["SWRabs_dir"]),
        FT(input_vars["SWRabs_diff"]),
        Bool(input_vars["RESPreCalc"]),
        rsRoofPreCalc,
    )

    @test results[1] ≈ output_vars["Hroof_imp"]
    @test results[2] ≈ output_vars["Hroof_veg"]
    @test results[3] ≈ output_vars["Eroof_imp"]
    @test results[4] ≈ output_vars["Eroof_veg"]
    @test results[5] ≈ output_vars["Eroof_ground"]
    @test results[6] ≈ output_vars["Eroof_soil"]
    @test results[7] ≈ output_vars["TEroof_veg"]
    @test results[8] ≈ output_vars["LEroof_imp"]
    @test results[9] ≈ output_vars["LEroof_veg"]
    @test results[10] ≈ output_vars["LEroof_ground"]
    @test results[11] ≈ output_vars["LEroof_soil"]
    @test results[12] ≈ output_vars["LTEroof_veg"]
    @test results[13] ≈ output_vars["Ci_sun_roof"]
    @test results[14] ≈ output_vars["Ci_shd_roof"]
    @test results[15] ≈ output_vars["ra"]
    @test results[16] ≈ output_vars["rb"]
    @test isnan(results[17])  # rap_L
    @test results[18] ≈ output_vars["r_soil"]
    @test results[19] ≈ output_vars["rs_sun"]
    @test results[20] ≈ output_vars["rs_shd"]
end
