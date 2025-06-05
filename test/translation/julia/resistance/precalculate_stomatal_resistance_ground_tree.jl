using Test
using MAT
using UrbanTethysChloris.Resistance: precalculate_stomatal_resistance_ground_tree
using UrbanTethysChloris.ModelComponents.Parameters:
    UrbanGeometryParameters,
    LocationSpecificSurfaceFractions,
    VegetatedOpticalProperties,
    SimpleOpticalProperties,
    HeightDependentVegetationParameters,
    WindowParameters
using UrbanTethysChloris.ModelComponents.ForcingInputs:
    SunPositionInputs, MeteorologicalInputs
using UrbanTethysChloris.RayTracing: ViewFactor
using ....TestUtils:
    create_urban_geometry_parameters,
    create_location_specific_surface_fractions,
    create_vegetated_optical_properties,
    create_simple_optical_properties,
    create_height_dependent_vegetation_parameters,
    create_window_parameters

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "resistance_functions.PrecalculateStomatalResistanceGroundTree.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

# Convert Dict{String, Any} to named tuples
TempVec_ittm = (;
    TCanyon=input_vars["TempVec_ittm"]["TCanyon"],
    TGroundVeg=input_vars["TempVec_ittm"]["TGroundVeg"],
    TTree=input_vars["TempVec_ittm"]["TTree"],
)

Humidity_ittm = (; CanyonSpecific=input_vars["Humidity_ittm"]["CanyonSpecific"])

ParVegGround = create_height_dependent_vegetation_parameters(
    FT;
    LAI=input_vars["ParVegGround"]["LAI"],
    Kopt=input_vars["ParVegGround"]["Kopt"],
    Knit=input_vars["ParVegGround"]["Knit"],
    Psi_sto_50=input_vars["ParVegGround"]["Psi_sto_50"],
    Psi_sto_00=input_vars["ParVegGround"]["Psi_sto_00"],
    CT=Int64(input_vars["ParVegGround"]["CT"]),
    Vmax=input_vars["ParVegGround"]["Vmax"],
    DSE=input_vars["ParVegGround"]["DSE"],
    Ha=input_vars["ParVegGround"]["Ha"],
    FI=input_vars["ParVegGround"]["FI"],
    Do=input_vars["ParVegGround"]["Do"],
    a1=input_vars["ParVegGround"]["a1"],
    go=input_vars["ParVegGround"]["go"],
    e_rel=input_vars["ParVegGround"]["e_rel"],
    e_relN=input_vars["ParVegGround"]["e_relN"],
    gmes=input_vars["ParVegGround"]["gmes"],
    rjv=input_vars["ParVegGround"]["rjv"],
    mSl=input_vars["ParVegGround"]["mSl"],
    Sl=input_vars["ParVegGround"]["Sl"],
)

SoilPotW_ittm = (;
    SoilPotWGroundVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWGroundVeg_L"],
    SoilPotWGroundTot_H=input_vars["SoilPotW_ittm"]["SoilPotWGroundTot_H"],
)

CiCO2Leaf_ittm = (;
    CiCO2LeafTreeSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeSun"],
    CiCO2LeafTreeShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeShd"],
    CiCO2LeafGroundVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegSun"],
    CiCO2LeafGroundVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegShd"],
)

MeteoData = (;
    Pre=input_vars["MeteoData"]["Pre"],
    Catm_O2=input_vars["MeteoData"]["Catm_O2"],
    Catm_CO2=input_vars["MeteoData"]["Catm_CO2"],
    SW_dir=input_vars["MeteoData"]["SW_dir"],
    SW_diff=input_vars["MeteoData"]["SW_diff"],
)

geometry = create_urban_geometry_parameters(
    FT;
    hcanyon=input_vars["geometry"]["hcanyon"],
    wcanyon=input_vars["geometry"]["wcanyon"],
    htree=input_vars["geometry"]["htree"],
    radius_tree=radius_tree=input_vars["geometry"]["radius_tree"],
    distance_tree=input_vars["geometry"]["distance_tree"],
    trees=Bool(input_vars["ParTree"]["trees"]),
    ftree=input_vars["ParTree"]["ftree"],
)

FractionsGround = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsGround"]["fveg"],
    fbare=input_vars["FractionsGround"]["fbare"],
    fimp=input_vars["FractionsGround"]["fimp"],
)

PropOpticalGround = create_vegetated_optical_properties(
    FT;
    aveg=input_vars["PropOpticalGround"]["aveg"],
    abare=input_vars["PropOpticalGround"]["abare"],
    aimp=input_vars["PropOpticalGround"]["aimp"],
)

PropOpticalWall = create_simple_optical_properties(
    FT; albedo=input_vars["PropOpticalWall"]["albedo"]
)
PropOpticalTree = create_simple_optical_properties(
    FT; albedo=input_vars["PropOpticalTree"]["albedo"]
)

ParVegTree = create_height_dependent_vegetation_parameters(
    FT;
    LAI=input_vars["ParVegTree"]["LAI"],
    Kopt=input_vars["ParVegTree"]["Kopt"],
    Knit=input_vars["ParVegTree"]["Knit"],
    Psi_sto_50=input_vars["ParVegTree"]["Psi_sto_50"],
    Psi_sto_00=input_vars["ParVegTree"]["Psi_sto_00"],
    CT=Int64(input_vars["ParVegTree"]["CT"]),
    Vmax=input_vars["ParVegTree"]["Vmax"],
    DSE=input_vars["ParVegTree"]["DSE"],
    Ha=input_vars["ParVegTree"]["Ha"],
    FI=input_vars["ParVegTree"]["FI"],
    Do=input_vars["ParVegTree"]["Do"],
    a1=input_vars["ParVegTree"]["a1"],
    go=input_vars["ParVegTree"]["go"],
    e_rel=input_vars["ParVegTree"]["e_rel"],
    e_relN=input_vars["ParVegTree"]["e_relN"],
    gmes=input_vars["ParVegTree"]["gmes"],
    rjv=input_vars["ParVegTree"]["rjv"],
    mSl=input_vars["ParVegTree"]["mSl"],
    Sl=input_vars["ParVegTree"]["Sl"],
)

SunPosition = (;
    theta_n=input_vars["SunPosition"]["theta_n"],
    theta_Z=input_vars["SunPosition"]["theta_Z"],
)

View_Factor = ViewFactor{FT}(;
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

@testset "Zurich" begin
    rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L = precalculate_stomatal_resistance_ground_tree(
        TempVec_ittm,
        Humidity_ittm,
        ParVegGround,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        MeteoData,
        geometry,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParVegTree,
        SunPosition,
        View_Factor,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        input_vars["rb_L"],
        input_vars["rb_H"],
        input_vars["rap_can"],
        input_vars["rap_Htree_In"],
    )

    @test rs_sun_H ≈ output_vars["rs_sun_H"]
    @test rs_shd_H ≈ output_vars["rs_shd_H"]
    @test Ci_sun_H ≈ output_vars["Ci_sun_H"]
    @test Ci_shd_H ≈ output_vars["Ci_shd_H"]
    @test rs_sun_L ≈ output_vars["rs_sun_L"]
    @test rs_shd_L ≈ output_vars["rs_shd_L"]
    @test Ci_sun_L ≈ output_vars["Ci_sun_L"]
    @test Ci_shd_L ≈ output_vars["Ci_shd_L"]
end
