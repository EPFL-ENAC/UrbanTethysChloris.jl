using Test
using MAT
using UrbanTethysChloris.TurbulentHeat: heat_flux_ground
using ....TestUtils:
    create_height_dependent_vegetation_parameters,
    create_location_specific_surface_fractions,
    create_urban_geometry_parameters,
    create_vegetated_soil_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("turbulent_heat_function.HeatFlux_ground.json")

# Create structured inputs with correct types
FractionsGround = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsGround"]["fveg"],
    fbare=FT(input_vars["FractionsGround"]["fbare"]),
    fimp=input_vars["FractionsGround"]["fimp"],
)

Gemeotry_m = create_urban_geometry_parameters(
    FT;
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
    Height_tree=input_vars["Gemeotry_m"]["Height_tree"],
    Radius_tree=input_vars["Gemeotry_m"]["Radius_tree"],
    Hcan_max=input_vars["Gemeotry_m"]["Hcan_max"],
    Hcan_std=input_vars["Gemeotry_m"]["Hcan_std"],
    wcanyon=FT(input_vars["geometry"]["wcanyon"]),
    wroof_norm=input_vars["geometry"]["wroof_norm"],
    radius_tree=input_vars["geometry"]["radius_tree"],
    trees=Bool(input_vars["ParTree"]["trees"]),
)

ParVegGround = create_height_dependent_vegetation_parameters(
    FT;
    LAI=FT(input_vars["ParVegGround"]["LAI"]),
    SAI=input_vars["ParVegGround"]["SAI"],
    hc=input_vars["ParVegGround"]["hc"],
    d_leaf=input_vars["ParVegGround"]["d_leaf"],
    Kopt=input_vars["ParVegGround"]["Kopt"],
    Knit=input_vars["ParVegGround"]["Knit"],
    Psi_sto_50=FT(input_vars["ParVegGround"]["Psi_sto_50"]),
    Psi_sto_00=input_vars["ParVegGround"]["Psi_sto_00"],
    CT=Int(input_vars["ParVegGround"]["CT"]),
    Vmax=FT(input_vars["ParVegGround"]["Vmax"]),
    DSE=FT(input_vars["ParVegGround"]["DSE"]),
    Ha=FT(input_vars["ParVegGround"]["Ha"]),
    FI=input_vars["ParVegGround"]["FI"],
    Do=FT(input_vars["ParVegGround"]["Do"]),
    a1=FT(input_vars["ParVegGround"]["a1"]),
    go=input_vars["ParVegGround"]["go"],
    e_rel=FT(input_vars["ParVegGround"]["e_rel"]),
    e_relN=FT(input_vars["ParVegGround"]["e_relN"]),
    gmes=input_vars["ParVegGround"]["gmes"],
    rjv=input_vars["ParVegGround"]["rjv"],
    mSl=FT(input_vars["ParVegGround"]["mSl"]),
    Sl=input_vars["ParVegGround"]["Sl"],
    CASE_ROOT=Int(input_vars["ParVegGround"]["CASE_ROOT"]),
    ZR95=FT(input_vars["ParVegGround"]["ZR95"]),
    ZR50=input_vars["ParVegGround"]["ZR50"],
    ZRmax=input_vars["ParVegGround"]["ZRmax"],
)

ParVegTree = create_height_dependent_vegetation_parameters(
    FT;
    LAI=FT(input_vars["ParVegTree"]["LAI"]),
    SAI=input_vars["ParVegTree"]["SAI"],
    d_leaf=FT(input_vars["ParVegTree"]["d_leaf"]),
    Kopt=input_vars["ParVegTree"]["Kopt"],
    Knit=input_vars["ParVegTree"]["Knit"],
    Psi_sto_50=input_vars["ParVegTree"]["Psi_sto_50"],
    Psi_sto_00=input_vars["ParVegTree"]["Psi_sto_00"],
    CT=Int(input_vars["ParVegTree"]["CT"]),
    Vmax=FT(input_vars["ParVegTree"]["Vmax"]),
    DSE=input_vars["ParVegTree"]["DSE"],
    Ha=FT(input_vars["ParVegTree"]["Ha"]),
    FI=input_vars["ParVegTree"]["FI"],
    Do=FT(input_vars["ParVegTree"]["Do"]),
    a1=FT(input_vars["ParVegTree"]["a1"]),
    go=input_vars["ParVegTree"]["go"],
    e_rel=FT(input_vars["ParVegTree"]["e_rel"]),
    e_relN=FT(input_vars["ParVegTree"]["e_relN"]),
    gmes=input_vars["ParVegTree"]["gmes"],
    rjv=input_vars["ParVegTree"]["rjv"],
    mSl=FT(input_vars["ParVegTree"]["mSl"]),
    Sl=input_vars["ParVegTree"]["Sl"],
    SPARTREE=Int(input_vars["ParVegTree"]["SPARTREE"]),
    CASE_ROOT=Int(input_vars["ParVegTree"]["CASE_ROOT"]),
    ZR95=FT(input_vars["ParVegTree"]["ZR95"]),
    ZR50=input_vars["ParVegTree"]["ZR50"],
    ZRmax=input_vars["ParVegTree"]["ZRmax"],
)

ParSoilGround = create_vegetated_soil_parameters(
    FT;
    Pcla=input_vars["ParSoilGround"]["Pcla"],
    Psan=input_vars["ParSoilGround"]["Psan"],
    Porg=input_vars["ParSoilGround"]["Porg"],
    Kfc=input_vars["ParSoilGround"]["Kfc"],
    Phy=FT(input_vars["ParSoilGround"]["Phy"]),
    SPAR=Int(input_vars["ParSoilGround"]["SPAR"]),
    Kbot=input_vars["ParSoilGround"]["Kbot"],
    Sp_In=input_vars["ParSoilGround"]["Sp_In"],
    Zs=FT.(vec(input_vars["ParSoilGround"]["Zs"])),
)

TempVec_ittm = (;
    TGroundVeg=input_vars["TempVec_ittm"]["TGroundVeg"],
    TTree=input_vars["TempVec_ittm"]["TTree"],
    TCanyon=input_vars["TempVec_ittm"]["TCanyon"],
)

MeteoData = (;
    Tatm=input_vars["MeteoData"]["Tatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
    Zatm=FT(input_vars["MeteoData"]["Zatm"]),
    Uatm=input_vars["MeteoData"]["Uatm"],
    Catm_O2=input_vars["MeteoData"]["Catm_O2"],
    Catm_CO2=input_vars["MeteoData"]["Catm_CO2"],
)

ParInterceptionTree = (; Sp_In=input_vars["ParInterceptionTree"]["Sp_In"])

SoilPotW_ittm = (;
    SoilPotWGroundVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWGroundVeg_L"],
    SoilPotWGroundTot_H=input_vars["SoilPotW_ittm"]["SoilPotWGroundTot_H"],
)

Owater_ittm = (;
    OwGroundSoilImp=vec(input_vars["Owater_ittm"]["OwGroundSoilImp"]),
    OwGroundSoilBare=vec(input_vars["Owater_ittm"]["OwGroundSoilBare"]),
    OwGroundSoilVeg=vec(input_vars["Owater_ittm"]["OwGroundSoilVeg"]),
)

Vwater_ittm = (;
    VGroundSoilImp=vec(input_vars["Vwater_ittm"]["VGroundSoilImp"]),
    VGroundSoilBare=vec(input_vars["Vwater_ittm"]["VGroundSoilBare"]),
    VGroundSoilVeg=vec(input_vars["Vwater_ittm"]["VGroundSoilVeg"]),
)

ExWater_ittm = (;
    ExWaterGroundImp_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundImp_H"]),
    ExWaterGroundBare_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundBare_H"]),
    ExWaterGroundVeg_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_H"]),
    ExWaterGroundVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_L"]),
)

Int_ittm = (;
    IntGroundImp=input_vars["Int_ittm"]["IntGroundImp"],
    IntGroundBare=input_vars["Int_ittm"]["IntGroundBare"],
    IntGroundVegPlant=input_vars["Int_ittm"]["IntGroundVegPlant"],
    IntGroundVegGround=input_vars["Int_ittm"]["IntGroundVegGround"],
    IntTree=input_vars["Int_ittm"]["IntTree"],
)

CiCO2Leaf_ittm = (;
    CiCO2LeafTreeSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeSun"],
    CiCO2LeafTreeShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeShd"],
    CiCO2LeafGroundVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegSun"],
    CiCO2LeafGroundVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegShd"],
)

rsGroundPreCalc = (;
    rs_sun_L=input_vars["rsGroundPreCalc"]["rs_sun_L"],
    rs_shd_L=input_vars["rsGroundPreCalc"]["rs_shd_L"],
    Ci_sun_L=input_vars["rsGroundPreCalc"]["Ci_sun_L"],
    Ci_shd_L=input_vars["rsGroundPreCalc"]["Ci_shd_L"],
)

rsTreePreCalc = (;
    rs_sun_H=input_vars["rsTreePreCalc"]["rs_sun_H"],
    rs_shd_H=input_vars["rsTreePreCalc"]["rs_shd_H"],
    Ci_sun_H=input_vars["rsTreePreCalc"]["Ci_sun_H"],
    Ci_shd_H=input_vars["rsTreePreCalc"]["Ci_shd_H"],
)

ParCalculation = (;
    dth=input_vars["ParCalculation"]["dth"], row=input_vars["ParCalculation"]["row"]
)

@testset "heat_flux_ground" begin
    results = heat_flux_ground(
        vec(input_vars["TemperatureC"]),
        TempVec_ittm,
        MeteoData,
        Gemeotry_m,
        FractionsGround,
        ParVegGround,
        ParVegTree,
        ParSoilGround,
        SoilPotW_ittm,
        Owater_ittm,
        Vwater_ittm,
        ExWater_ittm,
        Int_ittm,
        CiCO2Leaf_ittm,
        ParInterceptionTree,
        ParCalculation,
        FT(input_vars["SWRdir_abs_tree"]),
        FT(input_vars["SWRdiff_abs_tree"]),
        FT(input_vars["SWRdir_abs_groundveg"]),
        FT(input_vars["SWRdiff_abs_groundveg"]),
        Bool(input_vars["RESPreCalc"]),
        rsGroundPreCalc,
        rsTreePreCalc,
    )

    @test results.Himp ≈ output_vars["Himp"]
    @test results.Hbare ≈ output_vars["Hbare"]
    @test results.Hveg ≈ output_vars["Hveg"]
    @test results.Htree ≈ output_vars["Htree"]
    @test results.Eimp ≈ output_vars["Eimp"]
    @test results.Ebare_pond ≈ output_vars["Ebare_pond"]
    @test results.Ebare_soil ≈ output_vars["Ebare_soil"]
    @test results.Eveg_int ≈ output_vars["Eveg_int"]
    @test results.Eveg_pond ≈ output_vars["Eveg_pond"]
    @test results.Eveg_soil ≈ output_vars["Eveg_soil"]
    @test results.TEveg ≈ output_vars["TEveg"]
    @test results.Etree_int ≈ output_vars["Etree_int"]
    @test results.TEtree ≈ output_vars["TEtree"]
    @test results.Ebare ≈ output_vars["Ebare"]
    @test results.Eveg ≈ output_vars["Eveg"]
    @test results.Etree ≈ output_vars["Etree"]
    @test results.LEimp ≈ output_vars["LEimp"]
    @test results.LEbare_pond ≈ output_vars["LEbare_pond"]
    @test results.LEbare_soil ≈ output_vars["LEbare_soil"]
    @test results.LEveg_int ≈ output_vars["LEveg_int"]
    @test results.LEveg_pond ≈ output_vars["LEveg_pond"]
    @test results.LEveg_soil ≈ output_vars["LEveg_soil"]
    @test results.LTEveg ≈ output_vars["LTEveg"]
    @test results.LEtree_int ≈ output_vars["LEtree_int"]
    @test results.LTEtree ≈ output_vars["LTEtree"]
    @test results.LEbare ≈ output_vars["LEbare"]
    @test results.LEveg ≈ output_vars["LEveg"]
    @test results.LEtree ≈ output_vars["LEtree"]
    @test results.Ci_sun_H ≈ output_vars["Ci_sun_H"]
    @test results.Ci_shd_H ≈ output_vars["Ci_shd_H"]
    @test results.Ci_sun_L ≈ output_vars["Ci_sun_L"]
    @test results.Ci_shd_L ≈ output_vars["Ci_shd_L"]
    @test results.rap_can ≈ output_vars["rap_can"]
    @test results.rap_Htree_In ≈ output_vars["rap_Htree_In"]
    @test results.rb_H ≈ output_vars["rb_H"]
    @test results.rb_L ≈ output_vars["rb_L"]
    @test results.r_soil_bare ≈ output_vars["r_soil_bare"]
    @test results.r_soil_veg ≈ output_vars["r_soil_veg"]
    @test results.alp_soil_bare ≈ output_vars["alp_soil_bare"]
    @test results.alp_soil_veg ≈ output_vars["alp_soil_veg"]
    @test results.rs_sun_L ≈ output_vars["rs_sun_L"]
    @test results.rs_shd_L ≈ output_vars["rs_shd_L"]
    @test results.rs_sun_H ≈ output_vars["rs_sun_H"]
    @test results.rs_shd_H ≈ output_vars["rs_shd_H"]
    @test results.u_Hcan ≈ output_vars["u_Hcan"]
    @test results.u_Zref_und ≈ output_vars["u_Zref_und"]
    @test results.Fsun_L ≈ output_vars["Fsun_L"]
    @test results.Fshd_L ≈ output_vars["Fshd_L"]
    @test results.dw_L ≈ output_vars["dw_L"]
end
