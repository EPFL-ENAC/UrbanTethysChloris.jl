using Test
using MAT
using UrbanTethysChloris.Resistance: precalculate_for_faster_numerical_solution
using UrbanTethysChloris.ModelComponents.Parameters:
    UrbanGeometryParameters,
    LocationSpecificSurfaceFractions,
    VegetatedOpticalProperties,
    SimpleOpticalProperties,
    HeightDependentVegetationParameters,
    WindowParameters
using UrbanTethysChloris.ModelComponents.ForcingInputs: SunPositionInputs
using UrbanTethysChloris.RayTracing: ViewFactor

using ....TestUtils:
    create_urban_geometry_parameters,
    create_location_specific_surface_fractions,
    create_vegetated_optical_properties,
    create_simple_optical_properties,
    create_height_dependent_vegetation_parameters,
    create_window_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data(
    "resistance_functions.PrecalculateForFasterNumericalSolution.json"
)

# Convert Dict{String, Any} to named tuples
TempVec_ittm = (;
    TCanyon=input_vars["TempVec_ittm"]["TCanyon"],
    TGroundVeg=input_vars["TempVec_ittm"]["TGroundVeg"],
    TTree=input_vars["TempVec_ittm"]["TTree"],
    TRoofVeg=input_vars["TempVec_ittm"]["TRoofVeg"],
    Tatm=input_vars["TempVec_ittm"]["Tatm"],
)

Humidity_ittm = (; CanyonSpecific=input_vars["Humidity_ittm"]["CanyonSpecific"])

HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])

ParVegGround = create_height_dependent_vegetation_parameters(
    FT;
    LAI=FT(input_vars["ParVegGround"]["LAI"]),
    Kopt=input_vars["ParVegGround"]["Kopt"],
    Knit=input_vars["ParVegGround"]["Knit"],
    Psi_sto_50=FT(input_vars["ParVegGround"]["Psi_sto_50"]),
    Psi_sto_00=input_vars["ParVegGround"]["Psi_sto_00"],
    CT=Int64(input_vars["ParVegGround"]["CT"]),
    Vmax=FT(input_vars["ParVegGround"]["Vmax"]),
    DSE=input_vars["ParVegGround"]["DSE"],
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
)

ParVegRoof = create_height_dependent_vegetation_parameters(
    FT;
    LAI=FT(input_vars["ParVegRoof"]["LAI"]),
    Kopt=input_vars["ParVegRoof"]["Kopt"],
    Knit=input_vars["ParVegRoof"]["Knit"],
    Psi_sto_50=FT(input_vars["ParVegRoof"]["Psi_sto_50"]),
    Psi_sto_00=input_vars["ParVegRoof"]["Psi_sto_00"],
    CT=Int64(input_vars["ParVegRoof"]["CT"]),
    Vmax=FT(input_vars["ParVegRoof"]["Vmax"]),
    DSE=input_vars["ParVegRoof"]["DSE"],
    Ha=FT(input_vars["ParVegRoof"]["Ha"]),
    FI=input_vars["ParVegRoof"]["FI"],
    Do=FT(input_vars["ParVegRoof"]["Do"]),
    a1=FT(input_vars["ParVegRoof"]["a1"]),
    go=input_vars["ParVegRoof"]["go"],
    e_rel=FT(input_vars["ParVegRoof"]["e_rel"]),
    e_relN=FT(input_vars["ParVegRoof"]["e_relN"]),
    gmes=input_vars["ParVegRoof"]["gmes"],
    rjv=input_vars["ParVegRoof"]["rjv"],
    mSl=FT(input_vars["ParVegRoof"]["mSl"]),
    Sl=input_vars["ParVegRoof"]["Sl"],
)

SoilPotW_ittm = (;
    SoilPotWGroundVeg_L=FT(input_vars["SoilPotW_ittm"]["SoilPotWGroundVeg_L"]),
    SoilPotWGroundTot_H=FT(input_vars["SoilPotW_ittm"]["SoilPotWGroundTot_H"]),
    SoilPotWRoofVeg_L=FT(input_vars["SoilPotW_ittm"]["SoilPotWRoofVeg_L"]),
)

CiCO2Leaf_ittm = (;
    CiCO2LeafTreeSun=FT(input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeSun"]),
    CiCO2LeafTreeShd=FT(input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeShd"]),
    CiCO2LeafGroundVegSun=FT(input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegSun"]),
    CiCO2LeafGroundVegShd=FT(input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegShd"]),
    CiCO2LeafRoofVegSun=FT(input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegSun"]),
    CiCO2LeafRoofVegShd=FT(input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegShd"]),
)

MeteoData = (;
    Pre=input_vars["MeteoData"]["Pre"],
    Catm_O2=FT(input_vars["MeteoData"]["Catm_O2"]),
    Catm_CO2=FT(input_vars["MeteoData"]["Catm_CO2"]),
    SW_dir=FT(input_vars["MeteoData"]["SW_dir"]),
    SW_diff=FT(input_vars["MeteoData"]["SW_diff"]),
    Zatm=input_vars["MeteoData"]["Zatm"],
    Uatm=input_vars["MeteoData"]["Uatm"],
    ea=input_vars["MeteoData"]["ea"],
    Tatm=input_vars["MeteoData"]["Tatm"],
)

geometry = create_urban_geometry_parameters(
    FT;
    hcanyon=input_vars["geometry"]["hcanyon"],
    wcanyon=FT(input_vars["geometry"]["wcanyon"]),
    htree=input_vars["geometry"]["htree"],
    radius_tree=input_vars["geometry"]["radius_tree"],
    distance_tree=input_vars["geometry"]["distance_tree"],
    trees=Bool(input_vars["ParTree"]["trees"]),
    ftree=FT(input_vars["ParTree"]["ftree"]),
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
    Height_tree=input_vars["Gemeotry_m"]["Height_tree"],
    Radius_tree=input_vars["Gemeotry_m"]["Radius_tree"],
    Hcan_max=input_vars["Gemeotry_m"]["Hcan_max"],
    Hcan_std=input_vars["Gemeotry_m"]["Hcan_std"],
)

FractionsGround = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsGround"]["fveg"],
    fbare=FT(input_vars["FractionsGround"]["fbare"]),
    fimp=input_vars["FractionsGround"]["fimp"],
)

FractionsRoof = create_location_specific_surface_fractions(
    FT; fveg=input_vars["FractionsRoof"]["fveg"], fimp=input_vars["FractionsRoof"]["fimp"]
)

PropOpticalGround = create_vegetated_optical_properties(
    FT;
    aveg=input_vars["PropOpticalGround"]["aveg"],
    abare=input_vars["PropOpticalGround"]["abare"],
    aimp=input_vars["PropOpticalGround"]["aimp"],
)

PropOpticalRoof = create_vegetated_optical_properties(
    FT;
    aveg=input_vars["PropOpticalRoof"]["aveg"],
    aimp=input_vars["PropOpticalRoof"]["aimp"],
)

PropOpticalWall = create_simple_optical_properties(
    FT; albedo=input_vars["PropOpticalWall"]["albedo"]
)

PropOpticalTree = create_simple_optical_properties(
    FT; albedo=input_vars["PropOpticalTree"]["albedo"]
)

ParVegTree = create_height_dependent_vegetation_parameters(
    FT;
    LAI=FT(input_vars["ParVegTree"]["LAI"]),
    Kopt=input_vars["ParVegTree"]["Kopt"],
    Knit=input_vars["ParVegTree"]["Knit"],
    Psi_sto_50=FT(input_vars["ParVegTree"]["Psi_sto_50"]),
    Psi_sto_00=input_vars["ParVegTree"]["Psi_sto_00"],
    CT=Int64(input_vars["ParVegTree"]["CT"]),
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
)

SunPosition = (;
    theta_n=input_vars["SunPosition"]["theta_n"],
    theta_Z=input_vars["SunPosition"]["theta_Z"],
)

ViewFact = ViewFactor{FT}(;
    F_gs_nT=input_vars["ViewFactor"]["F_gs_nT"],
    F_gw_nT=input_vars["ViewFactor"]["F_gw_nT"],
    F_ww_nT=input_vars["ViewFactor"]["F_ww_nT"],
    F_wg_nT=input_vars["ViewFactor"]["F_wg_nT"],
    F_ws_nT=input_vars["ViewFactor"]["F_ws_nT"],
    F_sg_nT=input_vars["ViewFactor"]["F_sg_nT"],
    F_sw_nT=input_vars["ViewFactor"]["F_sw_nT"],
    F_gs_T=input_vars["ViewFactor"]["F_gs_T"],
    F_gt_T=input_vars["ViewFactor"]["F_gt_T"],
    F_gw_T=input_vars["ViewFactor"]["F_gw_T"],
    F_ww_T=input_vars["ViewFactor"]["F_ww_T"],
    F_wt_T=input_vars["ViewFactor"]["F_wt_T"],
    F_wg_T=input_vars["ViewFactor"]["F_wg_T"],
    F_ws_T=input_vars["ViewFactor"]["F_ws_T"],
    F_sg_T=input_vars["ViewFactor"]["F_sg_T"],
    F_sw_T=input_vars["ViewFactor"]["F_sw_T"],
    F_st_T=input_vars["ViewFactor"]["F_st_T"],
    F_tg_T=input_vars["ViewFactor"]["F_tg_T"],
    F_tw_T=input_vars["ViewFactor"]["F_tw_T"],
    F_ts_T=input_vars["ViewFactor"]["F_ts_T"],
    F_tt_T=input_vars["ViewFactor"]["F_tt_T"],
)

ParWindows = create_window_parameters(
    FT;
    GlazingRatio=input_vars["ParWindows"]["GlazingRatio"],
    SolarAlbedo=input_vars["ParWindows"]["SolarAlbedo"],
)

# Create RES named tuple with properly dimensioned arrays
RES = (;
    raCanyontoAtmOrig=zeros(FT, 2, 1, 1),
    raRooftoAtm=zeros(FT, 2, 1, 1),
    rb_LRoof=zeros(FT, 2, 1, 1),
    rb_LGround=zeros(FT, 2, 1, 1),
    rb_HGround=zeros(FT, 2, 1, 1),
    rap_can=zeros(FT, 2, 1, 1),
    rap_Htree_In=zeros(FT, 2, 1, 1),
)

@testset "First iteration" begin
    ittn = 1
    ittm = 1

    fconv, rsRoofPreCalc, rsGroundPreCalc, rsTreePreCalc = precalculate_for_faster_numerical_solution(
        ittn,
        ittm,
        TempVec_ittm,
        Humidity_ittm,
        ParVegGround,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        MeteoData,
        HumidityAtm,
        geometry,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParVegTree,
        SunPosition,
        ViewFact,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        ParVegRoof,
        PropOpticalRoof,
        FractionsRoof,
        RES,
    )

    @test fconv == 0.0
    @test rsRoofPreCalc.rs_sun ≈ output_vars["rsRoofPreCalc"]["rs_sun"]
    @test rsRoofPreCalc.rs_shd ≈ output_vars["rsRoofPreCalc"]["rs_shd"]
    @test rsRoofPreCalc.Ci_sun ≈ output_vars["rsRoofPreCalc"]["Ci_sun"]
    @test rsRoofPreCalc.Ci_shd ≈ output_vars["rsRoofPreCalc"]["Ci_shd"]

    @test rsGroundPreCalc.rs_sun_L ≈ output_vars["rsGroundPreCalc"]["rs_sun_L"]
    @test rsGroundPreCalc.rs_shd_L ≈ output_vars["rsGroundPreCalc"]["rs_shd_L"]
    @test rsGroundPreCalc.Ci_sun_L ≈ output_vars["rsGroundPreCalc"]["Ci_sun_L"]
    @test rsGroundPreCalc.Ci_shd_L ≈ output_vars["rsGroundPreCalc"]["Ci_shd_L"]

    @test rsTreePreCalc.rs_sun_H ≈ output_vars["rsTreePreCalc"]["rs_sun_H"]
    @test rsTreePreCalc.rs_shd_H ≈ output_vars["rsTreePreCalc"]["rs_shd_H"]
    @test rsTreePreCalc.Ci_sun_H ≈ output_vars["rsTreePreCalc"]["Ci_sun_H"]
    @test rsTreePreCalc.Ci_shd_H ≈ output_vars["rsTreePreCalc"]["Ci_shd_H"]
end
