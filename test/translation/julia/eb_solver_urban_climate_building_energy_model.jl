using Test
using MAT
using UrbanTethysChloris: eb_solver_urban_climate_building_energy_model
using UrbanTethysChloris.RayTracing: ViewFactor
using ...TestUtils:
    load_matlab_data,
    create_urban_geometry_parameters,
    create_location_specific_surface_fractions,
    create_vegetated_soil_parameters,
    create_vegetated_optical_properties,
    create_simple_optical_properties,
    create_location_specific_thermal_properties,
    create_height_dependent_vegetation_parameters,
    create_indoor_optical_properties,
    create_hvac_parameters,
    create_window_parameters,
    create_thermal_building

FT = Float64
input_vars, output_vars = load_matlab_data("EBSolver_UrbanClimateBuildingEnergyModel.json")

TempVec_ittm = (;
    TWallSun=input_vars["TempVec_ittm"]["TWallSun"],
    TWallShade=input_vars["TempVec_ittm"]["TWallShade"],
    TWallIntSun=input_vars["TempVec_ittm"]["TWallIntSun"],
    TWallIntShade=input_vars["TempVec_ittm"]["TWallIntShade"],
    TGroundImp=input_vars["TempVec_ittm"]["TGroundImp"],
    TGroundBare=input_vars["TempVec_ittm"]["TGroundBare"],
    TGroundVeg=input_vars["TempVec_ittm"]["TGroundVeg"],
    TTree=input_vars["TempVec_ittm"]["TTree"],
    TCanyon=input_vars["TempVec_ittm"]["TCanyon"],
    TRoofVeg=input_vars["TempVec_ittm"]["TRoofVeg"],
    TRoofIntVeg=input_vars["TempVec_ittm"]["TRoofIntVeg"],
    TRoofIntImp=input_vars["TempVec_ittm"]["TRoofIntImp"],
)

TempVecB_ittm = (;
    Tbin=input_vars["TempVecB_ittm"]["Tbin"],
    qbin=input_vars["TempVecB_ittm"]["qbin"],
    Tinground=input_vars["TempVecB_ittm"]["Tinground"],
    Tintmass=input_vars["TempVecB_ittm"]["Tintmass"],
    Tinwallsun=input_vars["TempVecB_ittm"]["Tinwallsun"],
    Tinwallshd=input_vars["TempVecB_ittm"]["Tinwallshd"],
    Tceiling=input_vars["TempVecB_ittm"]["Tceiling"],
)

Humidity_ittm = (; CanyonSpecific=input_vars["Humidity_ittm"]["CanyonSpecific"])

MeteoData = (;
    SW_dir=FT(input_vars["MeteoData"]["SW_dir"]),
    SW_diff=FT(input_vars["MeteoData"]["SW_diff"]),
    LWR=input_vars["MeteoData"]["LWR"],
    Rain=input_vars["MeteoData"]["Rain"],
    Tatm=input_vars["MeteoData"]["Tatm"],
    Pre=input_vars["MeteoData"]["Pre"],
    ea=input_vars["MeteoData"]["ea"],
    Zatm=FT(input_vars["MeteoData"]["Zatm"]),
    Uatm=input_vars["MeteoData"]["Uatm"],
    q_atm=input_vars["MeteoData"]["q_atm"],
    Catm_O2=FT(input_vars["MeteoData"]["Catm_O2"]),
    Catm_CO2=FT(input_vars["MeteoData"]["Catm_CO2"]),
)

Int_ittm = (;
    IntGroundImp=input_vars["Int_ittm"]["IntGroundImp"],
    IntGroundVegPlant=input_vars["Int_ittm"]["IntGroundVegPlant"],
    IntGroundVegGround=input_vars["Int_ittm"]["IntGroundVegGround"],
    IntTree=input_vars["Int_ittm"]["IntTree"],
    IntGroundBare=input_vars["Int_ittm"]["IntGroundBare"],
    IntRoofImp=input_vars["Int_ittm"]["IntRoofImp"],
    IntRoofVegPlant=input_vars["Int_ittm"]["IntRoofVegPlant"],
    IntRoofVegGround=input_vars["Int_ittm"]["IntRoofVegGround"],
)

ExWater_ittm = (;
    ExWaterGroundImp_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundImp_H"]),
    ExWaterGroundBare_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundBare_H"]),
    ExWaterGroundVeg_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_H"]),
    ExWaterGroundVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_L"]),
    ExWaterRoofVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterRoofVeg_L"]),
)

Vwater_ittm = (;
    VGroundSoilVeg=vec(input_vars["Vwater_ittm"]["VGroundSoilVeg"]),
    VGroundSoilImp=vec(input_vars["Vwater_ittm"]["VGroundSoilImp"]),
    VGroundSoilBare=vec(input_vars["Vwater_ittm"]["VGroundSoilBare"]),
    VRoofSoilVeg=vec(input_vars["Vwater_ittm"]["VRoofSoilVeg"]),
)

Owater_ittm = (;
    OwGroundSoilVeg=vec(input_vars["Owater_ittm"]["OwGroundSoilVeg"]),
    OwGroundSoilBare=vec(input_vars["Owater_ittm"]["OwGroundSoilBare"]),
    OwGroundSoilImp=vec(input_vars["Owater_ittm"]["OwGroundSoilImp"]),
    OwRoofSoilVeg=vec(input_vars["Owater_ittm"]["OwRoofSoilVeg"]),
)

SoilPotW_ittm = (;
    SoilPotWGroundVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWGroundVeg_L"],
    SoilPotWGroundTot_H=input_vars["SoilPotW_ittm"]["SoilPotWGroundTot_H"],
    SoilPotWRoofVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWRoofVeg_L"],
)

CiCO2Leaf_ittm = (;
    CiCO2LeafGroundVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegSun"],
    CiCO2LeafGroundVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegShd"],
    CiCO2LeafTreeSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeSun"],
    CiCO2LeafTreeShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeShd"],
    CiCO2LeafRoofVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegSun"],
    CiCO2LeafRoofVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafRoofVegShd"],
)

TempDamp_ittm = (;
    TDampGroundImp=input_vars["TempDamp_ittm"]["TDampGroundImp"],
    TDampGroundBare=input_vars["TempDamp_ittm"]["TDampGroundBare"],
    TDampGroundVeg=input_vars["TempDamp_ittm"]["TDampGroundVeg"],
    TDampGroundBuild=input_vars["TempDamp_ittm"]["TDampGroundBuild"],
)

view_factor = ViewFactor{FT}(;
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

# Create parameter structs from input data
Gemeotry_m = create_urban_geometry_parameters(
    FT;
    Height_canyon=input_vars["Gemeotry_m"]["Height_canyon"],
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
    Width_roof=input_vars["Gemeotry_m"]["Width_roof"],
    Height_tree=input_vars["Gemeotry_m"]["Height_tree"],
    Radius_tree=input_vars["Gemeotry_m"]["Radius_tree"],
    Hcan_max=input_vars["Gemeotry_m"]["Hcan_max"],
    Distance_tree=input_vars["Gemeotry_m"]["Distance_tree"],
    Hcan_std=input_vars["Gemeotry_m"]["Hcan_std"],
    hcanyon=FT(input_vars["geometry"]["hcanyon"]),
    wcanyon=FT(input_vars["geometry"]["wcanyon"]),
    wroof_norm=input_vars["geometry"]["wroof_norm"],
    radius_tree=input_vars["geometry"]["radius_tree"],
    trees=Bool(input_vars["ParTree"]["trees"]),
    distance_tree=input_vars["geometry"]["distance_tree"],
    wcanyon_norm=input_vars["geometry"]["wcanyon_norm"],
    ratio=input_vars["geometry"]["ratio"],
    ftree=FT(input_vars["ParTree"]["ftree"]),
)

FractionsGround = create_location_specific_surface_fractions(
    FT;
    fveg=FT(input_vars["FractionsGround"]["fveg"]),
    fbare=FT(input_vars["FractionsGround"]["fbare"]),
    fimp=FT(input_vars["FractionsGround"]["fimp"]),
)

FractionsRoof = create_location_specific_surface_fractions(
    FT;
    fveg=input_vars["FractionsRoof"]["fveg"],
    fimp=input_vars["FractionsRoof"]["fimp"],
    Per_runoff=FT(input_vars["FractionsRoof"]["Per_runoff"]),
)

WallLayers = (;
    dz1_wall=input_vars["WallLayers"]["dz1_wall"],
    dz2_wall=input_vars["WallLayers"]["dz2_wall"],
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

ParInterceptionTree = (; Sp_In=input_vars["ParInterceptionTree"]["Sp_In"],)

PropOpticalGround = create_vegetated_optical_properties(
    FT;
    aimp=input_vars["PropOpticalGround"]["aimp"],
    aveg=input_vars["PropOpticalGround"]["aveg"],
    abare=input_vars["PropOpticalGround"]["abare"],
    eimp=input_vars["PropOpticalGround"]["eimp"],
    eveg=input_vars["PropOpticalGround"]["eveg"],
    ebare=input_vars["PropOpticalGround"]["ebare"],
)

PropOpticalWall = create_simple_optical_properties(
    FT;
    emissivity=input_vars["PropOpticalWall"]["emissivity"],
    albedo=input_vars["PropOpticalWall"]["albedo"],
)

PropOpticalTree = create_simple_optical_properties(
    FT;
    emissivity=input_vars["PropOpticalTree"]["emissivity"],
    albedo=input_vars["PropOpticalTree"]["albedo"],
)

ParThermalGround = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalGround"]["cv_s_imp"],
    lan_dry=input_vars["ParThermalGround"]["lan_dry_imp"],
)

ParThermalWall = create_location_specific_thermal_properties(
    FT;
    cv_s=input_vars["ParThermalWall"]["cv_s"],
    lan_dry=input_vars["ParThermalWall"]["lan_dry"],
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
    ZR95=FT.([input_vars["ParVegGround"]["ZR95"]]),
    ZR50=FT.([input_vars["ParVegGround"]["ZR50"]]),
    ZRmax=FT.([input_vars["ParVegGround"]["ZRmax"]]),
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
    ZR95=FT.([input_vars["ParVegTree"]["ZR95"]]),
    ZR50=FT.([input_vars["ParVegTree"]["ZR50"]]),
    ZRmax=FT.([input_vars["ParVegTree"]["ZRmax"]]),
)

ParSoilRoof = create_vegetated_soil_parameters(
    FT;
    Pcla=input_vars["ParSoilRoof"]["Pcla"],
    Psan=input_vars["ParSoilRoof"]["Psan"],
    Porg=input_vars["ParSoilRoof"]["Porg"],
    Kfc=input_vars["ParSoilRoof"]["Kfc"],
    Phy=FT(input_vars["ParSoilRoof"]["Phy"]),
    SPAR=Int(input_vars["ParSoilRoof"]["SPAR"]),
    Kbot=input_vars["ParSoilRoof"]["Kbot"],
    Sp_In=input_vars["ParSoilRoof"]["Sp_In"],
    Zs=FT.(vec(input_vars["ParSoilRoof"]["Zs"])),
    dz1=input_vars["ParSoilRoof"]["dz1"],
    dz2=input_vars["ParSoilRoof"]["dz2"],
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
    LAI=FT(input_vars["ParVegRoof"]["LAI"]),
    SAI=FT(input_vars["ParVegRoof"]["SAI"]),
    hc=FT(input_vars["ParVegRoof"]["hc"]),
    d_leaf=FT(input_vars["ParVegRoof"]["d_leaf"]),
    Kopt=FT(input_vars["ParVegRoof"]["Kopt"]),
    Knit=FT(input_vars["ParVegRoof"]["Knit"]),
    Psi_sto_50=FT(input_vars["ParVegRoof"]["Psi_sto_50"]),
    Psi_sto_00=FT(input_vars["ParVegRoof"]["Psi_sto_00"]),
    CT=Int(input_vars["ParVegRoof"]["CT"]),
    Vmax=FT(input_vars["ParVegRoof"]["Vmax"]),
    DSE=FT(input_vars["ParVegRoof"]["DSE"]),
    Ha=FT(input_vars["ParVegRoof"]["Ha"]),
    FI=FT(input_vars["ParVegRoof"]["FI"]),
    Do=FT(input_vars["ParVegRoof"]["Do"]),
    a1=FT(input_vars["ParVegRoof"]["a1"]),
    go=FT(input_vars["ParVegRoof"]["go"]),
    e_rel=FT(input_vars["ParVegRoof"]["e_rel"]),
    e_relN=FT(input_vars["ParVegRoof"]["e_relN"]),
    gmes=FT(input_vars["ParVegRoof"]["gmes"]),
    rjv=FT(input_vars["ParVegRoof"]["rjv"]),
    mSl=FT(input_vars["ParVegRoof"]["mSl"]),
    Sl=FT(input_vars["ParVegRoof"]["Sl"]),
    Rrootl=FT.([input_vars["ParVegRoof"]["Rrootl"]]),
    PsiL50=FT.([input_vars["ParVegRoof"]["PsiL50"]]),
    PsiX50=FT.([input_vars["ParVegRoof"]["PsiX50"]]),
    CASE_ROOT=Int(input_vars["ParVegRoof"]["CASE_ROOT"]),
    ZR95=FT.([input_vars["ParVegRoof"]["ZR95"]]),
    ZR50=FT.([input_vars["ParVegRoof"]["ZR50"]]),
    ZRmax=FT.([input_vars["ParVegRoof"]["ZRmax"]]),
)

SunPosition = (;
    theta_n=input_vars["SunPosition"]["theta_n"],
    theta_Z=input_vars["SunPosition"]["theta_Z"],
)

HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])

Anthropogenic = (;
    Qf_canyon=input_vars["Anthropogenic"]["Qf_canyon"],
    Tb=input_vars["Anthropogenic"]["Tb"],
    Waterf_roof=input_vars["Anthropogenic"]["Waterf_roof"],
)

ParCalculation = (;
    dth=Int(input_vars["ParCalculation"]["dth"]),
    dts=Int(input_vars["ParCalculation"]["dts"]),
    row=FT(input_vars["ParCalculation"]["row"]),
)

PropOpticalIndoors = create_indoor_optical_properties(
    FT;
    abc=input_vars["PropOpticalIndoors"]["abc"],
    abw=input_vars["PropOpticalIndoors"]["abw"],
    abg=input_vars["PropOpticalIndoors"]["abg"],
    abm=input_vars["PropOpticalIndoors"]["abm"],
    ec=input_vars["PropOpticalIndoors"]["ec"],
    eg=input_vars["PropOpticalIndoors"]["eg"],
    ew=input_vars["PropOpticalIndoors"]["ew"],
    em=input_vars["PropOpticalIndoors"]["em"],
)

ParHVAC = create_hvac_parameters(
    FT;
    ACon=Bool(input_vars["ParHVAC"]["ACon"]),
    Heatingon=Bool(input_vars["ParHVAC"]["Heatingon"]),
    TsetpointCooling=input_vars["ParHVAC"]["TsetpointCooling"],
    TsetpointHeating=input_vars["ParHVAC"]["TsetpointHeating"],
    RHsetpointCooling=FT(input_vars["ParHVAC"]["RHsetpointCooling"]),
    COPAC=input_vars["ParHVAC"]["COPAC"],
    COPHeat=input_vars["ParHVAC"]["COPHeat"],
    ACH=input_vars["ParHVAC"]["ACH"],
)

ParThermalBulidFloor = create_thermal_building(
    FT;
    IntMassOn=Bool(input_vars["ParThermalBulidFloor"]["IntMassOn"]),
    lan_ground_floor=input_vars["ParThermalBulidFloor"]["lan_ground_floor"],
    cv_ground_floor=input_vars["ParThermalBulidFloor"]["cv_ground_floor"],
    cv_floor_IntMass=input_vars["ParThermalBulidFloor"]["cv_floor_IntMass"],
    cv_wall_IntMass=input_vars["ParThermalBulidFloor"]["cv_wall_IntMass"],
    dzFloor=input_vars["ParThermalBulidFloor"]["dzFloor"],
    dzWall=input_vars["ParThermalBulidFloor"]["dzWall"],
    FloorHeight=FT(input_vars["ParThermalBulidFloor"]["FloorHeight"]),
)

ParWindows = create_window_parameters(
    FT;
    WindowsOn=Int(input_vars["ParWindows"]["WindowsOn"]),
    GlazingRatio=input_vars["ParWindows"]["GlazingRatio"],
    Uvalue=input_vars["ParWindows"]["Uvalue"],
    cv_glass=input_vars["ParWindows"]["cv_glass"],
    dztot=input_vars["ParWindows"]["dztot"],
)

rsRoofPreCalc = (;
    rs_sun=input_vars["rsRoofPreCalc"]["rs_sun"],
    rs_shd=input_vars["rsRoofPreCalc"]["rs_shd"],
    Ci_sun=input_vars["rsRoofPreCalc"]["Ci_sun"],
    Ci_shd=input_vars["rsRoofPreCalc"]["Ci_shd"],
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

HVACSchedule = (;
    Hequip=FT(input_vars["HVACSchedule"]["Hequip"]),
    Hpeople=FT(input_vars["HVACSchedule"]["Hpeople"]),
    LEequip=FT(input_vars["HVACSchedule"]["LEequip"]),
    LEpeople=FT(input_vars["HVACSchedule"]["LEpeople"]),
    AirConRoomFraction=FT(input_vars["HVACSchedule"]["AirConRoomFraction"]),
)

@testset "MATLAB" begin
    Ytot = eb_solver_urban_climate_building_energy_model(
        vec(input_vars["TemperatureTot"]),
        TempVec_ittm,
        TempVecB_ittm,
        Humidity_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        TempDamp_ittm,
        view_factor,
        Gemeotry_m,
        FractionsGround,
        FractionsRoof,
        WallLayers,
        ParSoilGround,
        ParInterceptionTree,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParThermalGround,
        ParThermalWall,
        ParVegGround,
        ParVegTree,
        ParSoilRoof,
        PropOpticalRoof,
        ParThermalRoof,
        ParVegRoof,
        SunPosition,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        PropOpticalIndoors,
        ParHVAC,
        ParThermalBulidFloor,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        Bool(input_vars["RESPreCalc"]),
        FT(input_vars["fconvPreCalc"]),
        FT(input_vars["fconv"]),
        rsRoofPreCalc,
        rsGroundPreCalc,
        rsTreePreCalc,
        HVACSchedule,
    )

    @test Ytot ≈ vec(output_vars["Ytot"])
end
