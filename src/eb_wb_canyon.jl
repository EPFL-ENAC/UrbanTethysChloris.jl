"""
    eb_wb_canyon(
        TemperatureC::Vector{FT},
        TemperatureB::Vector{FT},
        TempVec_ittm::NamedTuple,
        Humidity_ittm::NamedTuple,
        MeteoData::NamedTuple,
        Int_ittm::NamedTuple,
        ExWater_ittm::NamedTuple,
        Vwater_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        SoilPotW_ittm::NamedTuple,
        CiCO2Leaf_ittm::NamedTuple,
        TempDamp_ittm::NamedTuple,
        ViewFactor::RayTracing.ViewFactor{FT},
        Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        WallLayers::NamedTuple,
        ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParInterceptionTree::NamedTuple,
        PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
        PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
        PropOpticalTree::ModelComponents.Parameters.SimpleOpticalProperties{FT},
        ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        SunPosition::NamedTuple,
        HumidityAtm::NamedTuple,
        Anthropogenic::NamedTuple,
        ParCalculation::NamedTuple,
        TempVecB_ittm::NamedTuple,
        G2Roof::FT,
        PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
        ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
        ParThermalBuildingFloor::ModelComponents.Parameters.ThermalBuilding{FT},
        ParWindows::ModelComponents.Parameters.WindowParameters{FT},
        BEM_on::Bool,
        RESPreCalc::Bool,
        fconvPreCalc::FT,
        fconv::FT,
        rsGroundPreCalc::NamedTuple,
        rsTreePreCalc::NamedTuple,
        HVACSchedule::NamedTuple,
    ) where {FT<:AbstractFloat}

Calculate energy balance for canyon surfaces.

# Arguments
- `TemperatureC`: Canyon temperatures vector
- `TemperatureB`: Building temperatures vector
- `TempVec_ittm`: Temperature vectors at previous time step
- `Humidity_ittm`: Humidity at previous time step
- `MeteoData`: Meteorological data
- `Int_ittm`: Previous timestep interception values
- `ExWater_ittm`: Previous timestep extractable water values
- `Vwater_ittm`: Previous timestep soil water volume values
- `Owater_ittm`: Previous timestep soil water content values
- `SoilPotW_ittm`: Previous timestep soil water potential values
- `CiCO2Leaf_ittm`: Previous timestep leaf CO2 concentration values
- `TempDamp_ittm`: Previous timestep ground dampening temperature
- `ViewFactor`: View factors between surfaces
- `Gemeotry_m`: Urban geometry parameters
- `FractionsGround`: Ground surface fractions
- `WallLayers`: Wall layer parameters
- `ParSoilGround`: Soil parameters for ground
- `ParInterceptionTree`: Tree interception parameters
- `PropOpticalGround`: Ground optical properties
- `PropOpticalWall`: Wall optical properties
- `PropOpticalTree`: Tree optical properties
- `ParThermalGround`: Ground thermal properties
- `ParThermalWall`: Wall thermal properties
- `ParVegGround`: Ground vegetation parameters
- `ParVegTree`: Tree vegetation parameters
- `SunPosition`: Solar position parameters
- `HumidityAtm`: Atmospheric humidity parameters
- `Anthropogenic`: Anthropogenic parameters
- `ParCalculation`: Calculation parameters
- `TempVecB_ittm`: Building temperature vectors at previous time step
- `G2Roof`: Conductive heat flux through roof [W/m²]
- `PropOpticalIndoors`: Indoor optical properties
- `ParHVAC`: HVAC parameters
- `ParThermalBuildingFloor`: Building floor thermal parameters
- `ParWindows`: Window parameters
- `BEM_on`: Building Energy Model switch
- `RESPreCalc`: Use pre-calculated resistances
- `fconvPreCalc`: Pre-calculated convection factors
- `fconv`: Convection factors
- `rsGroundPreCalc`: Pre-calculated ground resistances
- `rsTreePreCalc`: Pre-calculated tree resistances
- `HVACSchedule`: HVAC operation schedule

# Returns
- `Ycanyon::Vector{FT}`: Canyon energy balance residuals
- `G2WallSun::FT`: Conductive heat flux through sunlit wall
- `G2WallShade::FT`: Conductive heat flux through shaded wall
- `SWRabs_t::Radiation.RadiationFluxes{FT}`: Absorbed shortwave radiation
- `SWRabsWallSunTransmitted::FT`: Shortwave radiation absorbed by sunlit wall transmitted indoors
- `SWRabsWallShadeTransmitted::FT`: Shortwave radiation absorbed by shaded wall transmitted indoors
"""
function eb_wb_canyon!(
    model::Model{FT},
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    Humidity_ittm::ModelComponents.ModelVariables.Humidity{FT},
    ViewFactor::RayTracing.ViewFactor{FT},
    WallLayers::NamedTuple,
    ParInterceptionTree::NamedTuple,
    ParCalculation::NamedTuple,
    G2Roof::FT,
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    BEM_on::Bool,
    RESPreCalc::Bool,
    fconvPreCalc::Bool,
    fconv::FT,
    rsGroundPreCalc::NamedTuple,
    rsTreePreCalc::NamedTuple,
) where {FT<:AbstractFloat}
    results = eb_wb_canyon(
        TemperatureC,
        TemperatureB,
        TempVec_ittm,
        Humidity_ittm,
        model.forcing.meteorological,
        model.variables.waterflux.Interception,
        model.variables.waterflux.ExWater,
        model.variables.waterflux.Vwater,
        model.variables.waterflux.Owater,
        model.variables.waterflux.SoilPotW,
        model.variables.waterflux.CiCO2Leaf,
        model.variables.temperature.tempdamp,
        model.variables.waterflux.Runon,
        model.variables.waterflux.Qinlat,
        ViewFactor,
        model.parameters.urbangeometry,
        model.parameters.surfacefractions.ground,
        WallLayers,
        model.parameters.soil.ground,
        ParInterceptionTree,
        model.parameters.optical.ground,
        model.parameters.optical.wall,
        model.parameters.optical.tree,
        model.parameters.thermal.ground,
        model.parameters.thermal.wall,
        model.parameters.vegetation.ground,
        model.parameters.vegetation.tree,
        model.forcing.sunposition,
        model.forcing.meteorological,
        model.forcing.anthropogenic,
        ParCalculation,
        model.variables.buildingenergymodel.TempVecB,
        G2Roof,
        model.parameters.building_energy.indoor_optical,
        ParHVAC,
        model.parameters.building_energy.thermal,
        model.parameters.building_energy.windows,
        BEM_on,
        RESPreCalc,
        fconvPreCalc,
        fconv,
        rsGroundPreCalc,
        rsTreePreCalc,
        model.forcing.hvacschedule,
    )

    update!(model.variables, results, eb_wb_canyon_dispatcher)

    return results
end

function eb_wb_canyon(
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    Humidity_ittm::ModelComponents.ModelVariables.Humidity{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    Int_ittm::ModelComponents.ModelVariables.Interception{FT},
    ExWater_ittm::ModelComponents.ModelVariables.ExWater{FT,MR,MG},
    Vwater_ittm::ModelComponents.ModelVariables.Vwater{FT,MR,MG},
    Owater_ittm::ModelComponents.ModelVariables.Owater{FT,MR,MG},
    SoilPotW_ittm::ModelComponents.ModelVariables.SoilPotW{FT},
    CiCO2Leaf_ittm::ModelComponents.ModelVariables.CiCO2Leaf{FT},
    TempDamp_ittm::ModelComponents.ModelVariables.TempDamp{FT},
    Runon_ittm::ModelComponents.ModelVariables.Runon{FT},
    Qinlat_ittm::ModelComponents.ModelVariables.Qinlat{FT},
    ViewFactor::RayTracing.ViewFactor{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    WallLayers::NamedTuple,
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParInterceptionTree::NamedTuple,
    PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    PropOpticalTree::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SunPosition::ModelComponents.ForcingInputs.SunPositionInputs{FT},
    HumidityAtm::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    Anthropogenic::ModelComponents.ForcingInputs.AnthropogenicInputs{FT,0},
    ParCalculation::NamedTuple,
    TempVecB_ittm::ModelComponents.ModelVariables.TempVecB{FT},
    G2Roof::FT,
    PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    ParThermalBuildingFloor::ModelComponents.Parameters.ThermalBuilding{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
    RESPreCalc::Bool,
    fconvPreCalc::Bool,
    fconv::FT,
    rsGroundPreCalc::NamedTuple,
    rsTreePreCalc::NamedTuple,
    HVACSchedule::ModelComponents.ForcingInputs.HVACSchedule{FT,0};
) where {FT<:AbstractFloat,MR,MG}

    # Calculate shortwave radiation
    SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t, albedo_canyon = Radiation.total_shortwave_absorbed(
        Gemeotry_m,
        MeteoData.SW_dir,
        MeteoData.SW_diff,
        SunPosition.theta_n,
        SunPosition.theta_Z,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParVegTree,
        ViewFactor,
        ParWindows,
        BEM_on,
    )

    # Tree absorbed: conversion from sphere to horizontal projected area
    SWRabs_t = setproperties(
        SWRabs_t,
        (
            Tree=SWRabs_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )
    SWRabsDir_t = setproperties(
        SWRabsDir_t,
        (
            Tree=SWRabsDir_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )
    SWRabsDiff_t = setproperties(
        SWRabsDiff_t,
        (
            Tree=SWRabsDiff_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )

    # Handle window adjustments
    SWRabsWindowSun = zero(FT)
    SWRtransWindowSun = zero(FT)
    SWRabsWindowShade = zero(FT)
    SWRtransWindowShade = zero(FT)
    SWRabsWallShadeTransmitted = zero(FT)
    SWRabsWallSunTransmitted = zero(FT)
    if BEM_on
        GlazingRatio = ParWindows.WindowsOn == 1 ? ParWindows.GlazingRatio : zero(FT)

        SWRabsWallSunExt=(1-GlazingRatio) * SWRabs_t.WallSun
        SWRtransWindowSun=SWRabs_t.WallSun
        SWRabsWallSunTransmitted=GlazingRatio * SWRabs_t.WallSun
        SWRabsWallShadeExt=(1-GlazingRatio) * SWRabs_t.WallShade
        SWRtransWindowShade=SWRabs_t.WallShade
        SWRabsWallShadeTransmitted=GlazingRatio * SWRabs_t.WallShade

    else
        SWRabsWallSunExt=SWRabs_t.WallSun
        SWRabsWallShadeExt=SWRabs_t.WallShade
    end

    SWRabs_t2 = AbsorbedRadiationFluxes(
        SWRabs_t,
        SWRabsWindowSun,
        SWRtransWindowSun,
        SWRabsWindowShade,
        SWRtransWindowShade,
        SWRabsWallSunTransmitted,
        SWRabsWallShadeTransmitted,
        SWRabsWallSunExt,
        SWRabsWallShadeExt,
    )

    if abs(SWRabsWallSunExt+SWRabsWallSunTransmitted - SWRabs_t.WallSun) > 1e-8
        @warn "Warning"
    elseif abs(SWRabsWallShadeExt + SWRabsWallShadeTransmitted - SWRabs_t.WallShade) > 1e-8
        @warn "Warning"
    end

    # Calculate longwave radiation
    LWRin_t, LWRout_t, LWRabs_t, LWREB_t = Radiation.total_longwave_absorbed(
        TemperatureC,
        Gemeotry_m,
        MeteoData.LWR,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ViewFactor,
    )

    # Tree absorbed: conversion from sphere to horizontal projected area
    LWRabs_t = setproperties(
        LWRabs_t,
        (
            Tree=LWRabs_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )

    # Calculate conductive heat fluxes
    # Sunlit wall
    G1WallSun, G2WallSun, dsWallSun = ConductiveHeat.conductive_heat_flux_walls(
        TemperatureC,
        TemperatureB,
        TempVec_ittm,
        TempVecB_ittm,
        Anthropogenic,
        ParThermalWall,
        WallLayers,
        ParCalculation,
        true,
        ParWindows,
        BEM_on,
    )

    # Shaded wall
    G1WallShade, G2WallShade, dsWallShade = ConductiveHeat.conductive_heat_flux_walls(
        TemperatureC,
        TemperatureB,
        TempVec_ittm,
        TempVecB_ittm,
        Anthropogenic,
        ParThermalWall,
        WallLayers,
        ParCalculation,
        false,
        ParWindows,
        BEM_on,
    )

    # Ground conductive heat fluxes
    G1GroundImp, TDampGroundImp = ConductiveHeat.conductive_heat_flux_ground_fr(
        TemperatureC,
        TempDamp_ittm,
        TempVec_ittm,
        Owater_ittm,
        ParCalculation,
        ParThermalGround,
        FractionsGround,
        ParSoilGround,
        ParVegTree,
        ParVegGround,
    )

    # Bare ground
    G1GroundBare, TDampGroundBare = ConductiveHeat.conductive_heat_flux_ground_vb(
        TemperatureC,
        TempDamp_ittm,
        Owater_ittm,
        TempVec_ittm,
        ParCalculation,
        ParSoilGround,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        0,
    )

    # Vegetated ground
    G1GroundVeg, TDampGroundVeg = ConductiveHeat.conductive_heat_flux_ground_vb(
        TemperatureC,
        TempDamp_ittm,
        Owater_ittm,
        TempVec_ittm,
        ParCalculation,
        ParSoilGround,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        1,
    )

    G1Ground =
        G1GroundImp*FractionsGround.fimp +
        G1GroundBare*FractionsGround.fbare +
        G1GroundVeg*FractionsGround.fveg;
    GTree = FT(NaN);
    dsGroundImp = FT(NaN);
    dsGroundBare = FT(NaN);
    dsGroundVeg = FT(NaN);
    dsTree = FT(NaN);
    dsCanyonAir = FT(NaN);
    TDampTree = FT(NaN);
    G1Canyon =
        Gemeotry_m.wcanyon*G1Ground +
        Gemeotry_m.hcanyon*G1WallSun +
        Gemeotry_m.hcanyon*G1WallShade;
    G2Canyon =
        Gemeotry_m.wcanyon*0 + Gemeotry_m.hcanyon*G2WallSun + Gemeotry_m.hcanyon*G2WallShade;

    # Calculate sensible and latent heat fluxes
    HfluxCanyon, LEfluxCanyon, raCanyontoAtm, raCanyontoAtmOrig, fconv, HumidityCan = TurbulentHeat.heat_flux_canyon(
        TemperatureC, Gemeotry_m, MeteoData, ParVegTree, fconvPreCalc, fconv
    )

    # Ground and tree heat fluxes
    HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, EfluxGroundImp, EfluxGroundBarePond, EfluxGroundBareSoil, EfluxGroundVegInt, EfluxGroundVegPond, EfluxGroundVegSoil, TEfluxGroundVeg, EfluxTreeInt, TEfluxTree, Ebare, Eveg, Etree, LEfluxGroundImp, LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg, LEfluxTreeInt, LTEfluxTree, LEbare, LEveg, LEtree, CiCO2LeafTreeSun, CiCO2LeafTreeShd, CiCO2LeafGroundVegSun, CiCO2LeafGroundVegShd, rap_can, rap_Htree_In, rb_HGround, rb_LGround, r_soilGroundbare, r_soilGroundveg, alp_soil_bare, alp_soil_veg, rs_sunGround, rs_shdGround, rs_sunTree, rs_shdTree, u_Hcan, u_Zref_und, Fsun_L, Fshd_L, dw_L = TurbulentHeat.heat_flux_ground(
        TemperatureC,
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
        SWRabsDir_t.Tree,
        SWRabsDiff_t.Tree,
        SWRabsDir_t.GroundVeg,
        SWRabsDiff_t.GroundVeg,
        RESPreCalc,
        rsGroundPreCalc,
        rsTreePreCalc,
    )

    # Wall heat fluxes
    HfluxWallSun, HfluxWallShade, EfluxWallSun, EfluxWallShade, LEfluxWallSun, LEfluxWallShade, RES_w1, RES_w2, rap_W1_In, rap_W2_In, Hwsun1, Hwshade1, Hwsun2, Hwshade2, cp_atm, rho_atm, L_heat, Zp1, Zp2, rap_Zp1, rap_Zp2 = TurbulentHeat.heat_flux_wall(
        TemperatureC, Gemeotry_m, MeteoData, ParVegTree, ParVegGround, FractionsGround
    )

    HfluxGround =
        HfluxGroundImp*FractionsGround.fimp +
        HfluxGroundBare*FractionsGround.fbare +
        HfluxGroundVeg*FractionsGround.fveg;

    LEfluxGroundBare = LEfluxGroundBarePond + LEfluxGroundBareSoil;
    LEfluxGroundVeg =
        LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg;
    LEfluxGround =
        LEfluxGroundImp*FractionsGround.fimp +
        LEfluxGroundBare*FractionsGround.fbare +
        LEfluxGroundVeg*FractionsGround.fveg;
    LEfluxTree = LEfluxTreeInt + LTEfluxTree;

    EfluxGroundBare = EfluxGroundBarePond + EfluxGroundBareSoil;
    EfluxGroundVeg =
        EfluxGroundVegInt + EfluxGroundVegPond + EfluxGroundVegSoil + TEfluxGroundVeg;
    EfluxGround =
        EfluxGroundImp*FractionsGround.fimp +
        EfluxGroundBare*FractionsGround.fbare +
        EfluxGroundVeg*FractionsGround.fveg;
    EfluxTree = EfluxTreeInt + TEfluxTree;
    EfluxCanyon = FT(NaN);

    # Building energy model contribution
    if BEM_on
        SWRinWsun = SWRabsWallSunTransmitted
        SWRinWshd = SWRabsWallShadeTransmitted

        _, WasteHeat = BuildingEnergyModel.eb_solver_building(
            TemperatureC,
            TemperatureB,
            TempVecB_ittm,
            TempVec_ittm,
            Humidity_ittm,
            MeteoData,
            SWRinWsun,
            SWRinWshd,
            G2Roof,
            G2WallSun,
            G2WallShade,
            TempDamp_ittm,
            SWRabs_t,
            Gemeotry_m,
            PropOpticalIndoors,
            ParHVAC,
            ParCalculation,
            ParThermalBuildingFloor,
            ParWindows,
            BEM_on,
            HVACSchedule,
        )
    else
        WasteHeat = (;
            SensibleFromAC_Can=zero(FT),
            LatentFromAC_Can=zero(FT),
            WaterFromAC_Can=zero(FT),
            SensibleFromHeat_Can=zero(FT),
            LatentFromHeat_Can=zero(FT),
            SensibleFromVent_Can=zero(FT),
            LatentFromVent_Can=zero(FT),
            TotAnthInput_URB=zero(FT),
        )
    end

    # Canyon air heat storage
    Vcanyon = (Gemeotry_m.Width_canyon * Gemeotry_m.Height_canyon) / Gemeotry_m.Width_canyon
    # cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    # rho_atm = Pre / (287.04 * Tatm) * (1 - (ea/Pre) * (1 - 0.622))
    # L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    dS_H_air =
        Vcanyon * cp_atm * rho_atm * (TemperatureC[9] - TempVec_ittm.TCanyon) /
        ParCalculation.dts
    dS_LE_air =
        Vcanyon * rho_atm * L_heat * (TemperatureC[10] - Humidity_ittm.CanyonSpecific) /
        ParCalculation.dts

    # Initialize energy balance residuals
    Ycanyon = zeros(FT, 10)

    Cimp = FractionsGround.fimp > 0
    Cbare = FractionsGround.fbare > 0
    Cveg = FractionsGround.fveg > 0
    Ctree = Gemeotry_m.trees == 1

    # Ground energy balances
    if FractionsGround.fimp > 0
        Ycanyon[1] =
            SWRabs_t.GroundImp + LWRabs_t.GroundImp - G1GroundImp - HfluxGroundImp -
            LEfluxGroundImp
    else
        Ycanyon[1] = TemperatureC[1] - 273.15
    end

    # This is correct because we take the fbare== 0 route
    if FractionsGround.fbare > 0
        Ycanyon[2] =
            SWRabs_t.GroundBare + LWRabs_t.GroundBare - G1GroundBare - HfluxGroundBare -
            LEfluxGroundBarePond - LEfluxGroundBareSoil
    else
        Ycanyon[2] = TemperatureC[2] - 273.15
    end

    if FractionsGround.fveg > 0
        Ycanyon[3] =
            SWRabs_t.GroundVeg + LWRabs_t.GroundVeg - G1GroundVeg - HfluxGroundVeg -
            LEfluxGroundVegInt - LEfluxGroundVegPond - LEfluxGroundVegSoil -
            LTEfluxGroundVeg
    else
        Ycanyon[3] = TemperatureC[3] - 273.15
    end

    # Wall energy balances
    Ycanyon[4] = SWRabsWallSunExt + LWRabs_t.WallSun - G1WallSun - HfluxWallSun
    Ycanyon[5] = SWRabsWallShadeExt + LWRabs_t.WallShade - G1WallShade - HfluxWallShade

    # Tree energy balance
    if Gemeotry_m.trees > 0
        Ycanyon[6] = SWRabs_t.Tree + LWRabs_t.Tree - HfluxTree - LEfluxTreeInt - LTEfluxTree
    else
        Ycanyon[6] = TemperatureC[6] - 273.15
    end

    # Wall interior energy
    Ycanyon[7] = G1WallSun - G2WallSun - dsWallSun
    Ycanyon[8] = G1WallShade - G2WallShade - dsWallShade

    # Canyon air temperature and humidity energy balances
    Ycanyon[9] =
        Anthropogenic.Qf_canyon +
        Cimp * FractionsGround.fimp * HfluxGroundImp +
        Cbare * FractionsGround.fbare * HfluxGroundBare +
        Cveg * FractionsGround.fveg * HfluxGroundVeg +
        Gemeotry_m.hcanyon * (HfluxWallSun + HfluxWallShade) +
        Ctree * 4 * Gemeotry_m.radius_tree * HfluxTree - HfluxCanyon - dS_H_air +
        WasteHeat.SensibleFromVent_Can +
        WasteHeat.SensibleFromAC_Can +
        WasteHeat.SensibleFromHeat_Can

    Ycanyon[10] =
        Cimp * FractionsGround.fimp * LEfluxGroundImp +
        Cbare * FractionsGround.fbare * (LEfluxGroundBarePond + LEfluxGroundBareSoil) +
        Cveg *
        FractionsGround.fveg *
        (LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg) +
        Ctree * 4 * Gemeotry_m.radius_tree * (LEfluxTreeInt + LTEfluxTree) - LEfluxCanyon -
        dS_LE_air +
        WasteHeat.LatentFromVent_Can +
        WasteHeat.LatentFromAC_Can +
        WasteHeat.LatentFromHeat_Can

    EBGroundImp = Ycanyon[1]
    EBGroundBare = Ycanyon[2]
    EBGroundVeg = Ycanyon[3]
    EBTree = Ycanyon[6]
    EBWallSun = Ycanyon[4]
    EBWallShade = Ycanyon[5]
    EBWallSunInt = Ycanyon[7]
    EBWallShadeInt = Ycanyon[8]
    EBCanyonT = Ycanyon[9]
    EBCanyonQ = Ycanyon[10]

    # Air temperature at 2m
    T2m, DHi, Himp_2m, Hbare_2m, Hveg_2m, Hwsun_2m, Hwshade_2m, Hcan_2m = TurbulentHeat.calculate_t2m(
        TemperatureC[1],
        TemperatureC[2],
        TemperatureC[3],
        TemperatureC[4],
        TemperatureC[5],
        TemperatureC[9],
        Zp1,
        rap_Zp1,
        rap_W1_In,
        rb_LGround,
        RES_w1,
        FractionsGround,
        Gemeotry_m,
        ParVegGround,
        TempVec_ittm,
        cp_atm,
        rho_atm,
        ParCalculation,
        fconv,
        MeteoData,
    )

    # Air humidity at 2m
    function t2m_closure(q2m)
        return TurbulentHeat.air_humidity_2m(
            q2m,
            T2m,
            TemperatureC[1],
            TemperatureC[2],
            TemperatureC[3],
            TemperatureC[9],
            TemperatureC[10],
            rap_Zp1,
            rap_W1_In,
            rb_LGround,
            alp_soil_bare,
            r_soilGroundbare,
            alp_soil_veg,
            r_soilGroundveg,
            rs_sunGround,
            rs_shdGround,
            dw_L,
            Fsun_L,
            Fshd_L,
            FractionsGround,
            ParVegGround,
            EfluxGroundImp,
            Ebare,
            EfluxGroundVegInt,
            EfluxGroundVegPond,
            EfluxGroundVegSoil,
            TEfluxGroundVeg,
            MeteoData.Pre,
            Humidity_ittm,
            fconv,
            MeteoData,
            Gemeotry_m,
            rho_atm,
            Zp1,
            ParCalculation,
        )
    end

    q2m = find_zero(t2m_closure, TemperatureC[10]; xrtol=eps(FT))

    DEi, Eimp_2m, Ebare_soil_2m, Eveg_int_2m, Eveg_soil_2m, TEveg_2m, Ecan_2m, q2m, e_T2m, RH_T2m, qcan, e_Tcan, RH_Tcan = TurbulentHeat.air_humidity_2m_output(
        q2m,
        T2m,
        TemperatureC[1],
        TemperatureC[2],
        TemperatureC[3],
        TemperatureC[9],
        TemperatureC[10],
        rap_Zp1,
        rap_W1_In,
        rb_LGround,
        alp_soil_bare,
        r_soilGroundbare,
        alp_soil_veg,
        r_soilGroundveg,
        rs_sunGround,
        rs_shdGround,
        dw_L,
        Fsun_L,
        Fshd_L,
        FractionsGround,
        ParVegGround,
        EfluxGroundImp,
        Ebare,
        EfluxGroundVegInt,
        EfluxGroundVegPond,
        EfluxGroundVegSoil,
        TEfluxGroundVeg,
        MeteoData.Pre,
        Humidity_ittm,
        fconv,
        MeteoData,
        Gemeotry_m,
        rho_atm,
        Zp1,
        ParCalculation,
    )

    q_tree_dwn, In_tree, dIn_tree_dt, q_gveg_dwn, In_gveg, dIn_gveg_dt, q_gimp_runoff, In_gimp, dIn_gimp_dt, f_inf_gimp, q_gbare_runoff, In_gbare, dIn_gbare_dt, f_inf_gbare, q_gveg_runoff, In_gveg_pond, dIn_gveg_pond_dt, f_inf_gveg, V_gimp, O_gimp, OS_gimp, Lk_gimp, Psi_s_H_gimp, Psi_s_L_gimp, Exwat_H_gimp, Exwat_L_gimp, Rd_gimp, TEgveg_imp, TEtree_imp, Egimp_soil, dV_dt_gimp, Psi_Soil_gimp, Kf_gimp, V_gbare, O_gbare, OS_gbare, Lk_gbare, Psi_s_H_gbare, Psi_s_L_gbare, Exwat_H_gbare, Exwat_L_gbare, Rd_gbare, TEgveg_bare, TEtree_bare, Egbare_Soil, dV_dt_gbare, Psi_soil_gbare, Kf_gbare, V_gveg, O_gveg, OS_gveg, Lk_gveg, Psi_s_H_gveg, Psi_s_L_gveg, Exwat_H_gveg, Exwat_L_gveg, Rd_gveg, TEgveg_veg, TEtree_veg, Egveg_Soil, dV_dt_gveg, Psi_soil_gveg, Kf_gveg, Qin_imp, Qin_bare, Qin_veg, Qin_bare2imp, Qin_bare2veg, Qin_imp2bare, Qin_imp2veg, Qin_veg2imp, Qin_veg2bare, V, O, OS, Lk, Rd, dV_dt, Psi_s_L, Exwat_L, TEgveg_tot, Psi_s_H_tot, Exwat_H, TEtree_tot, EB_TEtree, EB_TEgveg, WBIndv, WBTot, Runoff, Runon_ittm, Etot, DeepGLk, StorageTot = Water.water_canyon(
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        Qinlat_ittm,
        EfluxTreeInt,
        EfluxGroundVegInt,
        EfluxGroundImp,
        EfluxGroundBarePond,
        EfluxGroundVegPond,
        EfluxGroundBareSoil,
        EfluxGroundVegSoil,
        TEfluxGroundVeg,
        TEfluxTree,
        ParSoilGround,
        ParInterceptionTree,
        ParCalculation,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        Gemeotry_m,
        Anthropogenic,
    )

    return (;
        SWRin_t,
        SWRout_t,
        SWRabs_t=SWRabs_t2,
        SWRabsDir_t,
        SWRabsDiff_t,
        SWREB_t,
        albedo_canyon,
        LWRin_t,
        LWRout_t,
        LWRabs_t,
        LWREB_t,
        HfluxGroundImp,
        HfluxGroundBare,
        HfluxGroundVeg,
        HfluxTree,
        HfluxGround,
        EfluxGroundImp,
        EfluxGroundBarePond,
        EfluxGroundBareSoil,
        EfluxGroundVegInt,
        EfluxGroundVegPond,
        EfluxGroundVegSoil,
        TEfluxGroundVeg,
        EfluxTreeInt,
        TEfluxTree,
        EfluxGroundBare,
        EfluxGroundVeg,
        EfluxGround,
        EfluxTree,
        LEfluxGroundImp,
        LEfluxGroundBarePond,
        LEfluxGroundBareSoil,
        LEfluxGroundVegInt,
        LEfluxGroundVegPond,
        LEfluxGroundVegSoil,
        LTEfluxGroundVeg,
        LEfluxTreeInt,
        LTEfluxTree,
        LEfluxGroundBare,
        LEfluxGroundVeg,
        LEfluxGround,
        LEfluxTree,
        CiCO2LeafTreeSun,
        CiCO2LeafTreeShd,
        CiCO2LeafGroundVegSun,
        CiCO2LeafGroundVegShd,
        raCanyontoAtm,
        raCanyontoAtmOrig,
        rap_can,
        rap_Htree_In,
        rb_HGround,
        rb_LGround,
        r_soilGroundbare,
        r_soilGroundveg,
        alp_soilGroundbare=alp_soil_bare,
        alp_soilGroundveg=alp_soil_veg,
        rs_sunGround,
        rs_shdGround,
        rs_sunTree,
        rs_shdTree,
        Fsun_L, # unused
        Fshd_L, # unused
        dw_L, # unused
        RES_w1,
        RES_w2,
        rap_W1_In,
        rap_W2_In,
        rap_Zp1,
        HfluxWallSun,
        HfluxWallShade,
        EfluxWallSun,
        EfluxWallShade,
        LEfluxWallSun,
        LEfluxWallShade,
        HfluxCanyon,
        LEfluxCanyon,
        EfluxCanyon,
        G1WallSun,
        G2WallSun,
        dsWallSun,
        G1WallShade,
        G2WallShade,
        dsWallShade,
        G1GroundImp,
        TDampGroundImp,
        G1GroundBare,
        TDampGroundBare,
        G1GroundVeg,
        TDampGroundVeg,
        GTree,
        TDampTree,
        G1Ground,
        G1Canyon,
        G2Canyon,
        dsGroundImp,
        dsGroundBare,
        dsGroundVeg,
        dsTree,
        dsCanyonAir,
        Ycanyon,
        QTree=q_tree_dwn,
        IntTree=In_tree,
        dInt_dtTree=dIn_tree_dt,
        QGroundVegDrip=q_gveg_dwn,
        IntGroundVegPlant=In_gveg,
        dInt_dtGroundVegPlant=dIn_gveg_dt,
        QGroundImp=q_gimp_runoff,
        IntGroundImp=In_gimp,
        dInt_dtGroundImp=dIn_gimp_dt,
        fGroundImp=f_inf_gimp,
        QGroundBarePond=q_gbare_runoff,
        IntGroundBare=In_gbare,
        dInt_dtGroundBare=dIn_gbare_dt,
        fGroundBare=f_inf_gbare,
        QGroundVegPond=q_gveg_runoff,
        IntGroundVegGround=In_gveg_pond,
        dInt_dtGroundVegGround=dIn_gveg_pond_dt,
        fGroundVeg=f_inf_gveg,
        VGroundSoilImp=V_gimp,
        OwGroundSoilImp=O_gimp,
        OSwGroundSoilImp=OS_gimp,
        LkGroundImp=Lk_gimp,
        SoilPotWGroundImp_H=Psi_s_H_gimp,
        SoilPotWGroundImp_L=Psi_s_L_gimp,
        ExWaterGroundImp_H=Exwat_H_gimp,
        ExWaterGroundImp_L=Exwat_L_gimp,
        Rd_gimp,
        TEgveg_imp,
        TEtree_imp,
        Egimp_soil,
        dVGroundSoilImp_dt=dV_dt_gimp,
        Psi_Soil_gimp,
        Kf_gimp,
        VGroundSoilBare=V_gbare,
        OwGroundSoilBare=O_gbare,
        OSwGroundSoilBare=OS_gbare,
        LkGroundBare=Lk_gbare,
        SoilPotWGroundBare_H=Psi_s_H_gbare,
        SoilPotWGroundBare_L=Psi_s_L_gbare,
        ExWaterGroundBare_H=Exwat_H_gbare,
        ExWaterGroundBare_L=Exwat_L_gbare,
        QGroundBareSoil=Rd_gbare,
        TEgveg_bare,
        TEtree_bare,
        Egbare_Soil,
        dVGroundSoilBare_dt=dV_dt_gbare,
        Psi_soil_gbare,
        Kf_gbare,
        VGroundSoilVeg=V_gveg,
        OwGroundSoilVeg=O_gveg,
        OSwGroundSoilVeg=OS_gveg,
        LkGroundVeg=Lk_gveg,
        SoilPotWGroundVeg_H=Psi_s_H_gveg,
        SoilPotWGroundVeg_L=Psi_s_L_gveg,
        ExWaterGroundVeg_H=Exwat_H_gveg,
        ExWaterGroundVeg_L=Exwat_L_gveg,
        QGroundVegSoil=Rd_gveg,
        TEgveg_veg,
        TEtree_veg,
        Egveg_Soil,
        dVGroundSoilVeg_dt=dV_dt_gveg,
        Psi_soil_gveg,
        Kf_gveg,
        Qin_imp,
        Qin_bare,
        Qin_veg,
        Qin_bare2imp,
        Qin_bare2veg,
        Qin_imp2bare,
        Qin_imp2veg,
        Qin_veg2imp,
        Qin_veg2bare,
        VGroundSoilTot=V,
        OwGroundSoilTot=O,
        OSwGroundSoilTot=OS,
        LkGround=Lk,
        Rd,
        dVGroundSoilTot_dt=dV_dt,
        SoilPotWGroundTot_L=Psi_s_L,
        ExWaterGroundTot_L=Exwat_L,
        TEgveg_tot,
        SoilPotWGroundTot_H=Psi_s_H_tot,
        ExWaterGroundTot_H=Exwat_H,
        TEtree_tot,
        EB_TEtree,
        EB_TEgveg,
        WBIndv, # NamedTuple
        WBTot, # NamedTuple
        RunoffGroundTot=Runoff,
        RunonGroundTot=Runon_ittm,
        Etot,
        DeepGLk,
        StorageTot,
        EBGroundImp,
        EBGroundBare,
        EBGroundVeg,
        EBTree,
        EBWallSun,
        EBWallShade,
        EBWallSunInt,
        EBWallShadeInt,
        EBCanyonT,
        EBCanyonQ,
        HumidityCan,
        HumidityAtm,
        u_Hcan,
        u_Zref_und,
        T2m,
        q2m,
        e_T2m,
        RH_T2m,
        qcan,
        e_Tcan,
        RH_Tcan,
        DHi,
        Himp_2m,
        Hbare_2m,
        Hveg_2m,
        Hwsun_2m,
        Hwshade_2m,
        Hcan_2m,
        DEi,
        Eimp_2m,
        Ebare_soil_2m,
        Eveg_int_2m,
        Eveg_soil_2m,
        TEveg_2m,
        Ecan_2m,
        dS_H_air,
        dS_LE_air,
    )
end

function eb_wb_canyon(
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::NamedTuple,
    Humidity_ittm::NamedTuple,
    MeteoData::NamedTuple,
    Int_ittm::NamedTuple,
    ExWater_ittm::NamedTuple,
    Vwater_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    SoilPotW_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    TempDamp_ittm::NamedTuple,
    Runon_ittm::NamedTuple,
    Qinlat_ittm::NamedTuple,
    ViewFactor::RayTracing.ViewFactor{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    WallLayers::NamedTuple,
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParInterceptionTree::NamedTuple,
    PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    PropOpticalTree::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SunPosition::NamedTuple,
    HumidityAtm::NamedTuple,
    Anthropogenic::NamedTuple,
    ParCalculation::NamedTuple,
    TempVecB_ittm::NamedTuple,
    G2Roof::FT,
    PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    ParThermalBuildingFloor::ModelComponents.Parameters.ThermalBuilding{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
    RESPreCalc::Bool,
    fconvPreCalc::FT,
    fconv::FT,
    rsGroundPreCalc::NamedTuple,
    rsTreePreCalc::NamedTuple,
    HVACSchedule::NamedTuple,
) where {FT<:AbstractFloat}

    # Calculate shortwave radiation
    SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t, albedo_canyon = Radiation.total_shortwave_absorbed(
        Gemeotry_m,
        MeteoData.SW_dir,
        MeteoData.SW_diff,
        SunPosition.theta_n,
        SunPosition.theta_Z,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParVegTree,
        ViewFactor,
        ParWindows,
        BEM_on,
    )

    # Tree absorbed: conversion from sphere to horizontal projected area
    SWRabs_t = setproperties(
        SWRabs_t,
        (
            Tree=SWRabs_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )
    SWRabsDir_t = setproperties(
        SWRabsDir_t,
        (
            Tree=SWRabsDir_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )
    SWRabsDiff_t = setproperties(
        SWRabsDiff_t,
        (
            Tree=SWRabsDiff_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )

    # Handle window adjustments
    SWRabsWindowSun = zero(FT)
    SWRtransWindowSun = zero(FT)
    SWRabsWindowShade = zero(FT)
    SWRtransWindowShade = zero(FT)
    SWRabsWallShadeTransmitted = zero(FT)
    SWRabsWallSunTransmitted = zero(FT)
    if BEM_on
        GlazingRatio = ParWindows.WindowsOn == 1 ? ParWindows.GlazingRatio : zero(FT)

        SWRabsWallSunExt=(1-GlazingRatio) * SWRabs_t.WallSun
        SWRtransWindowSun=SWRabs_t.WallSun
        SWRabsWallSunTransmitted=GlazingRatio * SWRabs_t.WallSun
        SWRabsWallShadeExt=(1-GlazingRatio) * SWRabs_t.WallShade
        SWRtransWindowShade=SWRabs_t.WallShade
        SWRabsWallShadeTransmitted=GlazingRatio * SWRabs_t.WallShade

    else
        SWRabsWallSunExt=SWRabs_t.WallSun
        SWRabsWallShadeExt=SWRabs_t.WallShade
    end

    if abs(SWRabsWallSunTransmitted - SWRabs_t.WallSun) > 1e-8
        @warn "Warning"
    elseif abs(SWRabsWallShadeExt + SWRabsWallShadeTransmitted - SWRabs_t.WallShade) > 1e-8
        @warn "Warning"
    end

    # Calculate longwave radiation
    LWRin_t, LWRout_t, LWRabs_t, LWREB_t = Radiation.total_longwave_absorbed(
        TemperatureC,
        Gemeotry_m,
        MeteoData.LWR,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ViewFactor,
    )

    # Tree absorbed: conversion from sphere to horizontal projected area
    LWRabs_t = setproperties(
        LWRabs_t,
        (
            Tree=LWRabs_t.Tree * 4 * Gemeotry_m.radius_tree * π /
                 (4 * Gemeotry_m.radius_tree),
        ),
    )

    # Calculate conductive heat fluxes
    # Sunlit wall
    G1WallSun, G2WallSun, dsWallSun = ConductiveHeat.conductive_heat_flux_walls(
        TemperatureC,
        TemperatureB,
        TempVec_ittm,
        TempVecB_ittm,
        Anthropogenic,
        ParThermalWall,
        WallLayers,
        ParCalculation,
        true,
        ParWindows,
        BEM_on,
    )

    # Shaded wall
    G1WallShade, G2WallShade, dsWallShade = ConductiveHeat.conductive_heat_flux_walls(
        TemperatureC,
        TemperatureB,
        TempVec_ittm,
        TempVecB_ittm,
        Anthropogenic,
        ParThermalWall,
        WallLayers,
        ParCalculation,
        false,
        ParWindows,
        BEM_on,
    )

    # Ground conductive heat fluxes
    G1GroundImp, TDampGroundImp = ConductiveHeat.conductive_heat_flux_ground_fr(
        TemperatureC,
        TempDamp_ittm,
        TempVec_ittm,
        Owater_ittm,
        ParCalculation,
        ParThermalGround,
        FractionsGround,
        ParSoilGround,
        ParVegTree,
        ParVegGround,
    )

    # Bare ground
    G1GroundBare, TDampGroundBare = ConductiveHeat.conductive_heat_flux_ground_vb(
        TemperatureC,
        TempDamp_ittm,
        Owater_ittm,
        TempVec_ittm,
        ParCalculation,
        ParSoilGround,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        0,
    )

    # Vegetated ground
    G1GroundVeg, TDampGroundVeg = ConductiveHeat.conductive_heat_flux_ground_vb(
        TemperatureC,
        TempDamp_ittm,
        Owater_ittm,
        TempVec_ittm,
        ParCalculation,
        ParSoilGround,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        1,
    )

    G1Ground =
        G1GroundImp*FractionsGround.fimp +
        G1GroundBare*FractionsGround.fbare +
        G1GroundVeg*FractionsGround.fveg;
    GTree = FT(NaN);
    dsGroundImp = FT(NaN);
    dsGroundBare = FT(NaN);
    dsGroundVeg = FT(NaN);
    dsTree = FT(NaN);
    dsCanyonAir = FT(NaN);
    TDampTree = FT(NaN);
    G1Canyon =
        Gemeotry_m.wcanyon*G1Ground +
        Gemeotry_m.hcanyon*G1WallSun +
        Gemeotry_m.hcanyon*G1WallShade;
    G2Canyon =
        Gemeotry_m.wcanyon*0 + Gemeotry_m.hcanyon*G2WallSun + Gemeotry_m.hcanyon*G2WallShade;

    # Calculate sensible and latent heat fluxes
    HfluxCanyon, LEfluxCanyon, raCanyontoAtm, raCanyontoAtmOrig, fconv, HumidityCan = TurbulentHeat.heat_flux_canyon(
        TemperatureC, Gemeotry_m, MeteoData, ParVegTree, fconvPreCalc, fconv
    )

    # Ground and tree heat fluxes
    HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, EfluxGroundImp, EfluxGroundBarePond, EfluxGroundBareSoil, EfluxGroundVegInt, EfluxGroundVegPond, EfluxGroundVegSoil, TEfluxGroundVeg, EfluxTreeInt, TEfluxTree, Ebare, Eveg, Etree, LEfluxGroundImp, LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg, LEfluxTreeInt, LTEfluxTree, LEbare, LEveg, LEtree, CiCO2LeafTreeSun, CiCO2LeafTreeShd, CiCO2LeafGroundVegSun, CiCO2LeafGroundVegShd, rap_can, rap_Htree_In, rb_HGround, rb_LGround, r_soilGroundbare, r_soilGroundveg, alp_soil_bare, alp_soil_veg, rs_sunGround, rs_shdGround, rs_sunTree, rs_shdTree, u_Hcan, u_Zref_und, Fsun_L, Fshd_L, dw_L = TurbulentHeat.heat_flux_ground(
        TemperatureC,
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
        SWRabsDir_t.Tree,
        SWRabsDiff_t.Tree,
        SWRabsDir_t.GroundVeg,
        SWRabsDiff_t.GroundVeg,
        RESPreCalc,
        rsGroundPreCalc,
        rsTreePreCalc,
    )

    # Wall heat fluxes
    HfluxWallSun, HfluxWallShade, EfluxWallSun, EfluxWallShade, LEfluxWallSun, LEfluxWallShade, RES_w1, RES_w2, rap_W1_In, rap_W2_In, Hwsun1, Hwshade1, Hwsun2, Hwshade2, cp_atm, rho_atm, L_heat, Zp1, Zp2, rap_Zp1, rap_Zp2 = TurbulentHeat.heat_flux_wall(
        TemperatureC, Gemeotry_m, MeteoData, ParVegTree, ParVegGround, FractionsGround
    )

    HfluxGround =
        HfluxGroundImp*FractionsGround.fimp +
        HfluxGroundBare*FractionsGround.fbare +
        HfluxGroundVeg*FractionsGround.fveg;

    LEfluxGroundBare = LEfluxGroundBarePond + LEfluxGroundBareSoil;
    LEfluxGroundVeg =
        LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg;
    LEfluxGround =
        LEfluxGroundImp*FractionsGround.fimp +
        LEfluxGroundBare*FractionsGround.fbare +
        LEfluxGroundVeg*FractionsGround.fveg;
    LEfluxTree = LEfluxTreeInt + LTEfluxTree;

    EfluxGroundBare = EfluxGroundBarePond + EfluxGroundBareSoil;
    EfluxGroundVeg =
        EfluxGroundVegInt + EfluxGroundVegPond + EfluxGroundVegSoil + TEfluxGroundVeg;
    EfluxGround =
        EfluxGroundImp*FractionsGround.fimp +
        EfluxGroundBare*FractionsGround.fbare +
        EfluxGroundVeg*FractionsGround.fveg;
    EfluxTree = EfluxTreeInt + TEfluxTree;
    EfluxCanyon = FT(NaN);

    # Building energy model contribution
    if BEM_on
        SWRinWsun = SWRabsWallSunTransmitted
        SWRinWshd = SWRabsWallShadeTransmitted

        _, WasteHeat = BuildingEnergyModel.eb_solver_building(
            TemperatureC,
            TemperatureB,
            TempVecB_ittm,
            TempVec_ittm,
            Humidity_ittm,
            MeteoData,
            SWRinWsun,
            SWRinWshd,
            G2Roof,
            G2WallSun,
            G2WallShade,
            TempDamp_ittm,
            SWRabs_t,
            Gemeotry_m,
            PropOpticalIndoors,
            ParHVAC,
            ParCalculation,
            ParThermalBuildingFloor,
            ParWindows,
            BEM_on,
            HVACSchedule,
        )
    else
        WasteHeat = (;
            SensibleFromAC_Can=zero(FT),
            LatentFromAC_Can=zero(FT),
            WaterFromAC_Can=zero(FT),
            SensibleFromHeat_Can=zero(FT),
            LatentFromHeat_Can=zero(FT),
            SensibleFromVent_Can=zero(FT),
            LatentFromVent_Can=zero(FT),
            TotAnthInput_URB=zero(FT),
        )
    end

    # Canyon air heat storage
    Vcanyon = (Gemeotry_m.Width_canyon * Gemeotry_m.Height_canyon) / Gemeotry_m.Width_canyon
    # cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    # rho_atm = Pre / (287.04 * Tatm) * (1 - (ea/Pre) * (1 - 0.622))
    # L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    dS_H_air =
        Vcanyon * cp_atm * rho_atm * (TemperatureC[9] - TempVec_ittm.TCanyon) /
        ParCalculation.dts
    dS_LE_air =
        Vcanyon * rho_atm * L_heat * (TemperatureC[10] - Humidity_ittm.CanyonSpecific) /
        ParCalculation.dts

    # Initialize energy balance residuals
    Ycanyon = zeros(FT, 10)

    Cimp = FractionsGround.fimp > 0
    Cbare = FractionsGround.fbare > 0
    Cveg = FractionsGround.fveg > 0
    Ctree = Gemeotry_m.trees == 1

    # Ground energy balances
    if FractionsGround.fimp > 0
        Ycanyon[1] =
            SWRabs_t.GroundImp + LWRabs_t.GroundImp - G1GroundImp - HfluxGroundImp -
            LEfluxGroundImp
    else
        Ycanyon[1] = TemperatureC[1] - 273.15
    end

    # This is correct because we take the fbare== 0 route
    if FractionsGround.fbare > 0
        Ycanyon[2] =
            SWRabs_t.GroundBare + LWRabs_t.GroundBare - G1GroundBare - HfluxGroundBare -
            LEfluxGroundBarePond - LEfluxGroundBareSoil
    else
        Ycanyon[2] = TemperatureC[2] - 273.15
    end

    if FractionsGround.fveg > 0
        Ycanyon[3] =
            SWRabs_t.GroundVeg + LWRabs_t.GroundVeg - G1GroundVeg - HfluxGroundVeg -
            LEfluxGroundVegInt - LEfluxGroundVegPond - LEfluxGroundVegSoil -
            LTEfluxGroundVeg
    else
        Ycanyon[3] = TemperatureC[3] - 273.15
    end

    # Wall energy balances
    Ycanyon[4] = SWRabs_t.WallSun + LWRabs_t.WallSun - G1WallSun - HfluxWallSun
    Ycanyon[5] = SWRabs_t.WallShade + LWRabs_t.WallShade - G1WallShade - HfluxWallShade

    # Tree energy balance
    if Gemeotry_m.trees > 0
        Ycanyon[6] = SWRabs_t.Tree + LWRabs_t.Tree - HfluxTree - LEfluxTreeInt - LTEfluxTree
    else
        Ycanyon[6] = TemperatureC[6] - 273.15
    end

    # Wall interior energy balances
    Ycanyon[7] = G1WallSun - G2WallSun - dsWallSun
    Ycanyon[8] = G1WallShade - G2WallShade - dsWallShade

    # Canyon air temperature and humidity energy balances
    Ycanyon[9] =
        Anthropogenic.Qf_canyon +
        Cimp * FractionsGround.fimp * HfluxGroundImp +
        Cbare * FractionsGround.fbare * HfluxGroundBare +
        Cveg * FractionsGround.fveg * HfluxGroundVeg +
        Gemeotry_m.hcanyon * (HfluxWallSun + HfluxWallShade) +
        Ctree * 4 * Gemeotry_m.radius_tree * HfluxTree - HfluxCanyon - dS_H_air +
        WasteHeat.SensibleFromVent_Can +
        WasteHeat.SensibleFromAC_Can +
        WasteHeat.SensibleFromHeat_Can

    Ycanyon[10] =
        Cimp * FractionsGround.fimp * LEfluxGroundImp +
        Cbare * FractionsGround.fbare * (LEfluxGroundBarePond + LEfluxGroundBareSoil) +
        Cveg *
        FractionsGround.fveg *
        (LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg) +
        Ctree * 4 * Gemeotry_m.radius_tree * (LEfluxTreeInt + LTEfluxTree) - LEfluxCanyon -
        dS_LE_air +
        WasteHeat.LatentFromVent_Can +
        WasteHeat.LatentFromAC_Can +
        WasteHeat.LatentFromHeat_Can

    EBGroundImp = Ycanyon[1]
    EBGroundBare = Ycanyon[2]
    EBGroundVeg = Ycanyon[3]
    EBTree = Ycanyon[6]
    EBWallSun = Ycanyon[4]
    EBWallShade = Ycanyon[5]
    EBWallSunInt = Ycanyon[7]
    EBWallShadeInt = Ycanyon[8]
    EBCanyonT = Ycanyon[9]
    EBCanyonQ = Ycanyon[10]

    # Air temperature at 2m
    T2m, DHi, Himp_2m, Hbare_2m, Hveg_2m, Hwsun_2m, Hwshade_2m, Hcan_2m = TurbulentHeat.calculate_t2m(
        TemperatureC[1],
        TemperatureC[2],
        TemperatureC[3],
        TemperatureC[4],
        TemperatureC[5],
        TemperatureC[9],
        Zp1,
        rap_Zp1,
        rap_W1_In,
        rb_LGround,
        RES_w1,
        FractionsGround,
        Gemeotry_m,
        ParVegGround,
        TempVec_ittm,
        cp_atm,
        rho_atm,
        ParCalculation,
        fconv,
        MeteoData,
    )

    # Air humidity at 2m
    function t2m_closure(q2m)
        return TurbulentHeat.air_humidity_2m(
            q2m,
            T2m,
            TemperatureC[1],
            TemperatureC[2],
            TemperatureC[3],
            TemperatureC[9],
            TemperatureC[10],
            rap_Zp1,
            rap_W1_In,
            rb_LGround,
            alp_soil_bare,
            r_soilGroundbare,
            alp_soil_veg,
            r_soilGroundveg,
            rs_sunGround,
            rs_shdGround,
            dw_L,
            Fsun_L,
            Fshd_L,
            FractionsGround,
            ParVegGround,
            EfluxGroundImp,
            Ebare,
            EfluxGroundVegInt,
            EfluxGroundVegPond,
            EfluxGroundVegSoil,
            TEfluxGroundVeg,
            MeteoData.Pre,
            Humidity_ittm,
            fconv,
            MeteoData,
            Gemeotry_m,
            rho_atm,
            Zp1,
            ParCalculation,
        )
    end

    q2m = find_zero(t2m_closure, TemperatureC[10]; xrtol=eps(FT))

    DEi, Eimp_2m, Ebare_soil_2m, Eveg_int_2m, Eveg_soil_2m, TEveg_2m, Ecan_2m, q2m, e_T2m, RH_T2m, qcan, e_Tcan, RH_Tcan = TurbulentHeat.air_humidity_2m_output(
        q2m,
        T2m,
        TemperatureC[1],
        TemperatureC[2],
        TemperatureC[3],
        TemperatureC[9],
        TemperatureC[10],
        rap_Zp1,
        rap_W1_In,
        rb_LGround,
        alp_soil_bare,
        r_soilGroundbare,
        alp_soil_veg,
        r_soilGroundveg,
        rs_sunGround,
        rs_shdGround,
        dw_L,
        Fsun_L,
        Fshd_L,
        FractionsGround,
        ParVegGround,
        EfluxGroundImp,
        Ebare,
        EfluxGroundVegInt,
        EfluxGroundVegPond,
        EfluxGroundVegSoil,
        TEfluxGroundVeg,
        MeteoData.Pre,
        Humidity_ittm,
        fconv,
        MeteoData,
        Gemeotry_m,
        rho_atm,
        Zp1,
        ParCalculation,
    )

    q_tree_dwn, In_tree, dIn_tree_dt, q_gveg_dwn, In_gveg, dIn_gveg_dt, q_gimp_runoff, In_gimp, dIn_gimp_dt, f_inf_gimp, q_gbare_runoff, In_gbare, dIn_gbare_dt, f_inf_gbare, q_gveg_runoff, In_gveg_pond, dIn_gveg_pond_dt, f_inf_gveg, V_gimp, O_gimp, OS_gimp, Lk_gimp, Psi_s_H_gimp, Psi_s_L_gimp, Exwat_H_gimp, Exwat_L_gimp, Rd_gimp, TEgveg_imp, TEtree_imp, Egimp_soil, dV_dt_gimp, Psi_Soil_gimp, Kf_gimp, V_gbare, O_gbare, OS_gbare, Lk_gbare, Psi_s_H_gbare, Psi_s_L_gbare, Exwat_H_gbare, Exwat_L_gbare, Rd_gbare, TEgveg_bare, TEtree_bare, Egbare_Soil, dV_dt_gbare, Psi_soil_gbare, Kf_gbare, V_gveg, O_gveg, OS_gveg, Lk_gveg, Psi_s_H_gveg, Psi_s_L_gveg, Exwat_H_gveg, Exwat_L_gveg, Rd_gveg, TEgveg_veg, TEtree_veg, Egveg_Soil, dV_dt_gveg, Psi_soil_gveg, Kf_gveg, Qin_imp, Qin_bare, Qin_veg, Qin_bare2imp, Qin_bare2veg, Qin_imp2bare, Qin_imp2veg, Qin_veg2imp, Qin_veg2bare, V, O, OS, Lk, Rd, dV_dt, Psi_s_L, Exwat_L, TEgveg_tot, Psi_s_H_tot, Exwat_H, TEtree_tot, EB_TEtree, EB_TEgveg, WBIndv, WBTot, Runoff, Runon_ittm, Etot, DeepGLk, StorageTot = Water.water_canyon(
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        Qinlat_ittm,
        EfluxTreeInt,
        EfluxGroundVegInt,
        EfluxGroundImp,
        EfluxGroundBarePond,
        EfluxGroundVegPond,
        EfluxGroundBareSoil,
        EfluxGroundVegSoil,
        TEfluxGroundVeg,
        TEfluxTree,
        ParSoilGround,
        ParInterceptionTree,
        ParCalculation,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        Gemeotry_m,
        Anthropogenic,
    )

    return (;
        SWRin_t,
        SWRout_t,
        SWRabs_t,
        SWRabsDir_t,
        SWRabsDiff_t,
        SWREB_t,
        albedo_canyon,
        LWRin_t,
        LWRout_t,
        LWRabs_t,
        LWREB_t,
        HfluxGroundImp,
        HfluxGroundBare,
        HfluxGroundVeg,
        HfluxTree,
        HfluxGround,
        EfluxGroundImp,
        EfluxGroundBarePond,
        EfluxGroundBareSoil,
        EfluxGroundVegInt,
        EfluxGroundVegPond,
        EfluxGroundVegSoil,
        TEfluxGroundVeg,
        EfluxTreeInt,
        TEfluxTree,
        EfluxGroundBare,
        EfluxGroundVeg,
        EfluxGround,
        EfluxTree,
        LEfluxGroundImp,
        LEfluxGroundBarePond,
        LEfluxGroundBareSoil,
        LEfluxGroundVegInt,
        LEfluxGroundVegPond,
        LEfluxGroundVegSoil,
        LTEfluxGroundVeg,
        LEfluxTreeInt,
        LTEfluxTree,
        LEfluxGroundBare,
        LEfluxGroundVeg,
        LEfluxGround,
        LEfluxTree,
        CiCO2LeafTreeSun,
        CiCO2LeafTreeShd,
        CiCO2LeafGroundVegSun,
        CiCO2LeafGroundVegShd,
        raCanyontoAtm,
        raCanyontoAtmOrig,
        rap_can,
        rap_Htree_In,
        rb_HGround,
        rb_LGround,
        r_soilGroundbare,
        r_soilGroundveg,
        alp_soil_bare,
        alp_soil_veg,
        rs_sunGround,
        rs_shdGround,
        rs_sunTree,
        rs_shdTree,
        Fsun_L,
        Fshd_L,
        dw_L,
        RES_w1,
        RES_w2,
        rap_W1_In,
        rap_W2_In,
        rap_Zp1,
        HfluxWallSun,
        HfluxWallShade,
        EfluxWallSun,
        EfluxWallShade,
        LEfluxWallSun,
        LEfluxWallShade,
        HfluxCanyon,
        LEfluxCanyon,
        EfluxCanyon,
        G1WallSun,
        G2WallSun,
        dsWallSun,
        G1WallShade,
        G2WallShade,
        dsWallShade,
        G1GroundImp,
        TDampGroundImp,
        G1GroundBare,
        TDampGroundBare,
        G1GroundVeg,
        TDampGroundVeg,
        GTree,
        TDampTree,
        G1Ground,
        G1Canyon,
        G2Canyon,
        dsGroundImp,
        dsGroundBare,
        dsGroundVeg,
        dsTree,
        dsCanyonAir,
        Ycanyon,
        q_tree_dwn,
        In_tree,
        dIn_tree_dt,
        q_gveg_dwn,
        In_gveg,
        dIn_gveg_dt,
        q_gimp_runoff,
        In_gimp,
        dIn_gimp_dt,
        f_inf_gimp,
        q_gbare_runoff,
        In_gbare,
        dIn_gbare_dt,
        f_inf_gbare,
        q_gveg_runoff,
        In_gveg_pond,
        dIn_gveg_pond_dt,
        f_inf_gveg,
        V_gimp,
        O_gimp,
        OS_gimp,
        Lk_gimp,
        Psi_s_H_gimp,
        Psi_s_L_gimp,
        Exwat_H_gimp,
        Exwat_L_gimp,
        Rd_gimp,
        TEgveg_imp,
        TEtree_imp,
        Egimp_soil,
        dV_dt_gimp,
        Psi_Soil_gimp,
        Kf_gimp,
        V_gbare,
        O_gbare,
        OS_gbare,
        Lk_gbare,
        Psi_s_H_gbare,
        Psi_s_L_gbare,
        Exwat_H_gbare,
        Exwat_L_gbare,
        Rd_gbare,
        TEgveg_bare,
        TEtree_bare,
        Egbare_Soil,
        dV_dt_gbare,
        Psi_soil_gbare,
        Kf_gbare,
        V_gveg,
        O_gveg,
        OS_gveg,
        Lk_gveg,
        Psi_s_H_gveg,
        Psi_s_L_gveg,
        Exwat_H_gveg,
        Exwat_L_gveg,
        Rd_gveg,
        TEgveg_veg,
        TEtree_veg,
        Egveg_Soil,
        dV_dt_gveg,
        Psi_soil_gveg,
        Kf_gveg,
        Qin_imp,
        Qin_bare,
        Qin_veg,
        Qin_bare2imp,
        Qin_bare2veg,
        Qin_imp2bare,
        Qin_imp2veg,
        Qin_veg2imp,
        Qin_veg2bare,
        V,
        O,
        OS,
        Lk,
        Rd,
        dV_dt,
        Psi_s_L,
        Exwat_L,
        TEgveg_tot,
        Psi_s_H_tot,
        Exwat_H,
        TEtree_tot,
        EB_TEtree,
        EB_TEgveg,
        WBIndv,
        WBTot,
        Runoff,
        Runon_ittm,
        Etot,
        DeepGLk,
        StorageTot,
        EBGroundImp,
        EBGroundBare,
        EBGroundVeg,
        EBTree,
        EBWallSun,
        EBWallShade,
        EBWallSunInt,
        EBWallShadeInt,
        EBCanyonT,
        EBCanyonQ,
        HumidityCan,
        HumidityAtm,
        u_Hcan,
        u_Zref_und,
        T2m,
        q2m,
        e_T2m,
        RH_T2m,
        qcan,
        e_Tcan,
        RH_Tcan,
        DHi,
        Himp_2m,
        Hbare_2m,
        Hveg_2m,
        Hwsun_2m,
        Hwshade_2m,
        Hcan_2m,
        DEi,
        Eimp_2m,
        Ebare_soil_2m,
        Eveg_int_2m,
        Eveg_soil_2m,
        TEveg_2m,
        Ecan_2m,
        dS_H_air,
        dS_LE_air,
    )
end
