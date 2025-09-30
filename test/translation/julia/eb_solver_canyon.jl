using Test
using MAT
using UrbanTethysChloris: eb_solver_canyon
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
input_vars, output_vars = load_matlab_data("EBSolver_canyon.json")

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

# Create additional parameter structs
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

ParInterceptionTree = (; Sp_In=input_vars["ParInterceptionTree"]["Sp_In"],)

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

ParWindows = create_window_parameters(
    FT;
    WindowsOn=Int(input_vars["ParWindows"]["WindowsOn"]),
    GlazingRatio=input_vars["ParWindows"]["GlazingRatio"],
    Uvalue=input_vars["ParWindows"]["Uvalue"],
    cv_glass=input_vars["ParWindows"]["cv_glass"],
    dztot=input_vars["ParWindows"]["dztot"],
)

# Create named tuples for other parameters
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
)

ExWater_ittm = (;
    ExWaterGroundImp_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundImp_H"]),
    ExWaterGroundBare_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundBare_H"]),
    ExWaterGroundVeg_H=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_H"]),
    ExWaterGroundVeg_L=vec(input_vars["ExWater_ittm"]["ExWaterGroundVeg_L"]),
)

Vwater_ittm = (;
    VGroundSoilVeg=vec(input_vars["Vwater_ittm"]["VGroundSoilVeg"]),
    VGroundSoilImp=vec(input_vars["Vwater_ittm"]["VGroundSoilImp"]),
    VGroundSoilBare=vec(input_vars["Vwater_ittm"]["VGroundSoilBare"]),
)

Owater_ittm = (;
    OwGroundSoilVeg=vec(input_vars["Owater_ittm"]["OwGroundSoilVeg"]),
    OwGroundSoilBare=vec(input_vars["Owater_ittm"]["OwGroundSoilBare"]),
    OwGroundSoilImp=vec(input_vars["Owater_ittm"]["OwGroundSoilImp"]),
)

SoilPotW_ittm = (;
    SoilPotWGroundVeg_L=input_vars["SoilPotW_ittm"]["SoilPotWGroundVeg_L"],
    SoilPotWGroundTot_H=input_vars["SoilPotW_ittm"]["SoilPotWGroundTot_H"],
)

CiCO2Leaf_ittm = (;
    CiCO2LeafGroundVegSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegSun"],
    CiCO2LeafGroundVegShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafGroundVegShd"],
    CiCO2LeafTreeSun=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeSun"],
    CiCO2LeafTreeShd=input_vars["CiCO2Leaf_ittm"]["CiCO2LeafTreeShd"],
)

TempDamp_ittm = (;
    TDampGroundImp=input_vars["TempDamp_ittm"]["TDampGroundImp"],
    TDampGroundBare=input_vars["TempDamp_ittm"]["TDampGroundBare"],
    TDampGroundVeg=input_vars["TempDamp_ittm"]["TDampGroundVeg"],
    TDampGroundBuild=input_vars["TempDamp_ittm"]["TDampGroundBuild"],
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

HumidityAtm = (; AtmVapourPreSat=input_vars["HumidityAtm"]["AtmVapourPreSat"])

Anthropogenic = (; Qf_canyon=input_vars["Anthropogenic"]["Qf_canyon"])

ParCalculation = (;
    dth=Int(input_vars["ParCalculation"]["dth"]),
    dts=Int(input_vars["ParCalculation"]["dts"]),
    row=FT(input_vars["ParCalculation"]["row"]),
)

HVACSchedule = (;
    Hequip=FT(input_vars["HVACSchedule"]["Hequip"]),
    Hpeople=FT(input_vars["HVACSchedule"]["Hpeople"]),
    LEequip=FT(input_vars["HVACSchedule"]["LEequip"]),
    LEpeople=FT(input_vars["HVACSchedule"]["LEpeople"]),
    AirConRoomFraction=FT(input_vars["HVACSchedule"]["AirConRoomFraction"]),
)

fconvPreCalc = input_vars["fconvPreCalc"]

fconv = input_vars["fconv"]

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

WallLayers = (;
    dz1_wall=input_vars["WallLayers"]["dz1_wall"],
    dz2_wall=input_vars["WallLayers"]["dz2_wall"],
)

SunPosition = (;
    theta_n=input_vars["SunPosition"]["theta_n"],
    theta_Z=input_vars["SunPosition"]["theta_Z"],
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

@testset "MATLAB" begin
    Ycanyon, G2WallSun, G2WallShade, SWRabs_t, SWRabsWallSunTransmitted, SWRabsWallShadeTransmitted = eb_solver_canyon(
        vec(input_vars["TemperatureC"]),
        vec(input_vars["TemperatureB"]),
        TempVec_ittm,
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
        SunPosition,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        TempVecB_ittm,
        FT(input_vars["G2Roof"]),
        PropOpticalIndoors,
        ParHVAC,
        ParThermalBulidFloor,
        ParWindows,
        Bool(input_vars["BEM_on"]),
        Bool(input_vars["RESPreCalc"]),
        FT(fconvPreCalc),
        FT(fconv),
        rsGroundPreCalc,
        rsTreePreCalc,
        HVACSchedule,
    )

    @test Ycanyon ≈ vec(output_vars["Ycanyon"])
    @test G2WallSun ≈ output_vars["G2WallSun"]
    @test G2WallShade ≈ output_vars["G2WallShade"]

    # Test shortwave radiation components
    @test SWRabs_t.GroundImp ≈ output_vars["SWRabs_t"]["SWRabsGroundImp"]
    @test SWRabs_t.GroundBare ≈ output_vars["SWRabs_t"]["SWRabsGroundBare"]
    @test SWRabs_t.GroundVeg ≈ output_vars["SWRabs_t"]["SWRabsGroundVeg"]
    @test SWRabs_t.WallSun ≈ output_vars["SWRabs_t"]["SWRabsWallSun"]
    @test SWRabs_t.WallShade ≈ output_vars["SWRabs_t"]["SWRabsWallShade"]
    @test SWRabs_t.Tree ≈ output_vars["SWRabs_t"]["SWRabsTree"]
    @test SWRabs_t.TotalGround ≈ output_vars["SWRabs_t"]["SWRabsTotalGround"]
    @test SWRabs_t.TotalCanyon ≈ output_vars["SWRabs_t"]["SWRabsTotalCanyon"]
end
