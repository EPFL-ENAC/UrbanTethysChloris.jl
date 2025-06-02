using Test
using MAT
using UrbanTethysChloris.Resistance: precalculate_stomatal_resistance_roof

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "resistance_functions.PrecalculateStomatalResistanceRoof.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

# convert Dict{String, Any} to a named tuple
TempVec_ittm = (; TRoofVeg=input_vars["TempVec_ittm"]["TRoofVeg"])
MeteoData = (;
    Tatm=input_vars["MeteoData"]["Tatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
    Catm_O2=input_vars["MeteoData"]["Catm_O2"],
    Catm_CO2=input_vars["MeteoData"]["Catm_CO2"],
    SW_dir=input_vars["MeteoData"]["SW_dir"],
    SW_diff=input_vars["MeteoData"]["SW_diff"],
)

HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])

ParVegRoof = (;
    LAI=input_vars["ParVegRoof"]["LAI"],
    Kopt=input_vars["ParVegRoof"]["Kopt"],
    Knit=input_vars["ParVegRoof"]["Knit"],
    Psi_sto_50=input_vars["ParVegRoof"]["Psi_sto_50"],
    Psi_sto_00=input_vars["ParVegRoof"]["Psi_sto_00"],
    CT=Int64(input_vars["ParVegRoof"]["CT"]),
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
)

SoilPotW_ittm = (; SoilPotWRoofVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWRoofVeg_L"])

CiCO2Leaf_ittm = (;
    CiCO2LeafRoofVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegSun"],
    CiCO2LeafRoofVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegShd"],
)

PropOpticalRoof = (; aveg=input_vars["PropOpticalRoof"]["aveg"],)

@testset "Zurich" begin
    rs_sun, rs_shd, Ci_sun, Ci_shd = precalculate_stomatal_resistance_roof(
        TempVec_ittm,
        MeteoData,
        HumidityAtm,
        ParVegRoof,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        PropOpticalRoof,
        input_vars["ra"],
        input_vars["rb"],
    )

    @test rs_sun ≈ output_vars["rs_sun"]
    @test rs_shd ≈ output_vars["rs_shd"]
    @test Ci_sun ≈ output_vars["Ci_sun"]
    @test Ci_shd ≈ output_vars["Ci_shd"]
end
