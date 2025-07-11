using Test
using MAT
using UrbanTethysChloris: eb_solver_roof
using ...TestUtils:
    create_height_dependent_vegetation_parameters,
    create_location_specific_surface_fractions,
    create_vegetated_soil_parameters,
    create_vegetated_optical_properties,
    create_location_specific_thermal_properties,
    create_urban_geometry_parameters

FT = Float64
dir = joinpath(@__DIR__, "..", "matlab", "data")
filename = "EBSolver_roof.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

# Create structured inputs with correct types
FractionsRoof = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsRoof"]["fveg"],
    fimp=input_vars["FractionsRoof"]["fimp"],
    Per_runoff=input_vars["FractionsRoof"]["Per_runoff"],
)

PropOpticalRoof = create_vegetated_optical_properties(
    FT;
    aimp=input_vars["PropOpticalRoof"]["aimp"],
    aveg=input_vars["PropOpticalRoof"]["aveg"],
    eimp=input_vars["PropOpticalRoof"]["eimp"],
    eveg=input_vars["PropOpticalRoof"]["eveg"],
)

ParThermalRoof = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalRoof"]["cv_s_imp"],
    lan_dry=input_vars["ParThermalRoof"]["lan_dry_imp"],
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
    Rrootl=[input_vars["ParVegRoof"]["Rrootl"]],
    PsiL50=[input_vars["ParVegRoof"]["PsiL50"]],
    PsiX50=[input_vars["ParVegRoof"]["PsiX50"]],
    CASE_ROOT=Int(input_vars["ParVegRoof"]["CASE_ROOT"]),
    ZR95=[input_vars["ParVegRoof"]["ZR95"]],
    ZR50=[input_vars["ParVegRoof"]["ZR50"]],
    ZRmax=[input_vars["ParVegRoof"]["ZRmax"]],
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
    dz1=input_vars["ParSoilRoof"]["dz1"],
    dz2=input_vars["ParSoilRoof"]["dz2"],
)

Gemeotry_m = create_urban_geometry_parameters(
    FT; Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"]
)

MeteoData = (;
    SW_dir=input_vars["MeteoData"]["SW_dir"],
    SW_diff=input_vars["MeteoData"]["SW_diff"],
    LWR=input_vars["MeteoData"]["LWR"],
    Rain=input_vars["MeteoData"]["Rain"],
    Tatm=input_vars["MeteoData"]["Tatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
    Zatm=input_vars["MeteoData"]["Zatm"],
    Uatm=input_vars["MeteoData"]["Uatm"],
    Catm_O2=input_vars["MeteoData"]["Catm_O2"],
    Catm_CO2=input_vars["MeteoData"]["Catm_CO2"],
)

TempVec_ittm = (;
    TRoofVeg=input_vars["TempVec_ittm"]["TRoofVeg"],
    TRoofIntVeg=input_vars["TempVec_ittm"]["TRoofIntVeg"],
    TRoofIntImp=input_vars["TempVec_ittm"]["TRoofIntImp"],
)
HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])
Anthropogenic = (;
    Tb=input_vars["Anthropogenic"]["Tb"],
    Waterf_roof=input_vars["Anthropogenic"]["Waterf_roof"],
)
ParCalculation = (;
    dth=Int(input_vars["ParCalculation"]["dth"]),
    dts=Int(input_vars["ParCalculation"]["dts"]),
    row=input_vars["ParCalculation"]["row"],
)

Int_ittm = (;
    IntRoofImp=input_vars["Int_ittm"]["IntRoofImp"],
    IntRoofVegPlant=input_vars["Int_ittm"]["IntRoofVegPlant"],
    IntRoofVegGround=input_vars["Int_ittm"]["IntRoofVegGround"],
)
ExWater_ittm = (; ExWaterRoofVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterRoofVeg_L"]))
Vwater_ittm = (; VRoofSoilVeg=vec(input_vars["Vwater_ittm"]["VRoofSoilVeg"]))
Owater_ittm = (; OwRoofSoilVeg=vec(input_vars["Owater_ittm"]["OwRoofSoilVeg"]))
SoilPotW_ittm = (; SoilPotWRoofVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWRoofVeg_L"])
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

@testset "MATLAB" begin
    Yroof, G2Roof = eb_solver_roof(
        vec(input_vars["TemperatureR"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        Gemeotry_m,
        FractionsRoof,
        ParSoilRoof,
        PropOpticalRoof,
        ParThermalRoof,
        ParVegRoof,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        Bool(input_vars["BEM_on"]),
        Bool(input_vars["RESPreCalc"]),
        rsRoofPreCalc,
    )

    @test Yroof ≈ vec(output_vars["Yroof"])
    @test G2Roof ≈ output_vars["G2Roof"]
end
