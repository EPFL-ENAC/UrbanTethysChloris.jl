"""
    eb_solver_canyon(
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
        ViewFactor::ViewFactor,
        Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        WallLayers::NamedTuple,
        ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParInterceptionTree::NamedTuple,
        PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
        PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
        PropOpticalTree::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
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
        fconvPreCalc::NamedTuple,
        fconv::NamedTuple,
        rsGroundPreCalc::NamedTuple,
        rsTreePreCalc::NamedTuple,
        HVACSchedule::NamedTuple
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
- `Gemoetry_m`: Urban geometry parameters
- `Gemoetry_m`: Canyon geometry parameters
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
- `Ycanyon`: Canyon energy balance residuals
- `G2WallSun`: Conductive heat flux through sunlit wall
- `G2WallShade`: Conductive heat flux through shaded wall
- `SWRabs_t`: Absorbed shortwave radiation
"""
function eb_solver_canyon(
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
    Gemoetry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    WallLayers::NamedTuple,
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParInterceptionTree::NamedTuple,
    PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    PropOpticalTree::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
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
    SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t = Radiation.total_shortwave_absorbed(
        Gemoetry_m,
        MeteoData.SWR_dir,
        MeteoData.SWR_diff,
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
    SWRabs_t = merge(
        SWRabs_t,
        (
            SWRabsTree=SWRabs_t.SWRabsTree * 4 * Gemoetry_m.radius_tree * π /
                       (4 * Gemoetry_m.radius_tree),
        ),
    )
    SWRabsDir_t = merge(
        SWRabsDir_t,
        (
            SWRabsTree=SWRabsDir_t.SWRabsTree * 4 * Gemoetry_m.radius_tree * π /
                       (4 * Gemoetry_m.radius_tree),
        ),
    )
    SWRabsDiff_t = merge(
        SWRabsDiff_t,
        (
            SWRabsTree=SWRabsDiff_t.SWRabsTree * 4 * Gemoetry_m.radius_tree * π /
                       (4 * Gemoetry_m.radius_tree),
        ),
    )

    # Handle window adjustments
    if BEM_on
        GlazingRatio = ParWindows.WindowsOn ? ParWindows.GlazingRatio : zero(FT)

        SWRabs_t = merge(
            SWRabs_t,
            (
                SWRabsWallSunExt=(1-GlazingRatio) * SWRabs_t.SWRabsWallSun,
                SWRabsWindowSun=zero(FT),
                SWRtransWindowSun=SWRabs_t.SWRabsWallSun,
                SWRabsWallSunTransmitted=GlazingRatio * SWRabs_t.SWRabsWallSun,
                SWRabsWallShadeExt=(1-GlazingRatio) * SWRabs_t.SWRabsWallShade,
                SWRabsWindowShade=zero(FT),
                SWRtransWindowShade=SWRabs_t.SWRabsWallShade,
                SWRabsWallShadeTransmitted=GlazingRatio * SWRabs_t.SWRabsWallShade,
            ),
        )
    else
        SWRabs_t = merge(
            SWRabs_t,
            (
                SWRabsWallSunExt=SWRabs_t.SWRabsWallSun,
                SWRabsWindowSun=zero(FT),
                SWRtransWindowSun=zero(FT),
                SWRabsWallSunTransmitted=zero(FT),
                SWRabsWallShadeExt=SWRabs_t.SWRabsWallShade,
                SWRabsWindowShade=zero(FT),
                SWRtransWindowShade=zero(FT),
                SWRabsWallShadeTransmitted=zero(FT),
            ),
        )
    end

    # Calculate longwave radiation
    LWRin_t, LWRout_t, LWRabs_t, LWREB_t = Radiation.total_lwr_absorbed(
        TemperatureC,
        Gemoetry_m,
        MeteoData,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ViewFactor,
    )

    # Tree absorbed: conversion from sphere to horizontal projected area
    LWRabs_t = merge(
        LWRabs_t,
        (
            LWRabsTree=LWRabs_t.LWRabsTree * 4 * Gemoetry_m.radius_tree * π /
                       (4 * Gemoetry_m.radius_tree),
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
    G1GroundImp, Tdp_ground_imp = ConductiveHeat.conductive_heat_flux_ground_fr(
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
    G1GroundBare, Tdp_ground_bare = ConductiveHeat.conductive_heat_flux_ground_vb(
        TemperatureC,
        TempDamp_ittm,
        Owater_ittm,
        TempVec_ittm,
        ParCalculation,
        ParSoilGround,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        false,
    )

    # Vegetated ground
    G1GroundVeg, Tdp_ground_veg = conductive_heat_flux_ground_vb(
        TemperatureC,
        TempDamp_ittm,
        Owater_ittm,
        TempVec_ittm,
        ParCalculation,
        ParSoilGround,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        true,
    )

    # Calculate sensible and latent heat fluxes
    HfluxCanyon, LEfluxCanyon, ra_canyon, ra_orig, fconv, HumidityCan = TurbulentHeat.heat_flux_canyon(
        TemperatureC, Gemoetry_m, MeteoData, ParVegTree, fconvPreCalc, fconv
    )

    # Ground and tree heat fluxes
    HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, Eground_imp_pond, Eground_bare_pond, Eground_bare_soil, Eground_veg_int, Eground_veg_pond, Eground_veg_soil, TEground_veg, E_tree_int, TE_tree, Ebare, Eveg, Etree, LEfluxGroundImp, LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg, LEfluxTreeInt, LTEfluxTree, LEbare, LEveg, LEtree, Ci_sun_tree, Ci_shd_tree, Ci_sun_ground, Ci_shd_ground, rap_can, rap_Htree_In, rb_H, rb_L, r_soil_bare, r_soil_veg, alp_soil_bare, alp_soil_veg, rs_sun_L, rs_shd_L, rs_sun_H, rs_shd_H, u_Hcan, u_Zref_und, Fsun_L, Fshd_L, dw_L = TurbulentHeat.heat_flux_ground(
        TemperatureC,
        TempVec_ittm,
        MeteoData,
        Gemoetry_m,
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
        SWRabsDir_t.SWRabsTree,
        SWRabsDiff_t.SWRabsTree,
        SWRabsDir_t.SWRabsGroundVeg,
        SWRabsDiff_t.SWRabsGroundVeg,
        RESPreCalc,
        rsGroundPreCalc,
        rsTreePreCalc,
    )

    # Wall heat fluxes
    HfluxWallSun, HfluxWallShade, Ewsun, Ewshade, LEwsun, LEwshade = TurbulentHeat.heat_flux_wall(
        TemperatureC, Gemoetry_m, MeteoData, ParVegTree, ParVegGround, FractionsGround
    )

    # Building energy model contribution
    if BEM_on
        SWRinWsun = SWRabs_t.SWRabsWallSunTransmitted
        SWRinWshd = SWRabs_t.SWRabsWallShadeTransmitted

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
            Gemoetry_m,
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
    Vcanyon = (Gemoetry_m.Width_canyon * Gemoetry_m.Height_canyon) / Gemoetry_m.Width_canyon
    cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    rho_atm = Pre / (287.04 * Tatm) * (1 - (ea/Pre) * (1 - 0.622))
    L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    dS_H_air =
        Vcanyon * cp_atm * rho_atm * (TemperatureC[9] - TempVec_ittm.TCanyon) /
        ParCalculation.dts
    dS_LE_air =
        Vcanyon * rho_atm * L_heat * (TemperatureC[10] - Humidity_ittm.CanyonSpecific) /
        ParCalculation.dts

    # Initialize energy balance residuals
    Ycanyon = zeros(FT, 10)

    # Energy balance calculations
    if Gemeotry_m.trees == 0
        SWRabs_t = merge(SWRabs_t, (SWRabsTree=zero(FT),))
        LWRabs_t = merge(LWRabs_t, (LWRabsTree=zero(FT),))
    end

    Cimp = FractionsGround.fimp > 0
    Cbare = FractionsGround.fbare > 0
    Cveg = FractionsGround.fveg > 0
    Ctree = Gemeotry_m.trees == 1

    # Ground energy balances
    if Cimp
        Ycanyon[1] =
            SWRabs_t.SWRabsGroundImp + LWRabs_t.LWRabsGroundImp - G1GroundImp -
            HfluxGroundImp - LEfluxGroundImp
    else
        Ycanyon[1] = TemperatureC[1] - 273.15
    end

    if Cbare
        Ycanyon[2] =
            SWRabs_t.SWRabsGroundBare + LWRabs_t.LWRabsGroundBare - G1GroundBare -
            HfluxGroundBare - LEfluxGroundBarePond - LEfluxGroundBareSoil
    else
        Ycanyon[2] = TemperatureC[2] - 273.15
    end

    if Cveg
        Ycanyon[3] =
            SWRabs_t.SWRabsGroundVeg + LWRabs_t.LWRabsGroundVeg - G1GroundVeg -
            HfluxGroundVeg - LEfluxGroundVegInt - LEfluxGroundVegPond -
            LEfluxGroundVegSoil - LTEfluxGroundVeg
    else
        Ycanyon[3] = TemperatureC[3] - 273.15
    end

    # Wall energy balances
    Ycanyon[4] =
        SWRabs_t.SWRabsWallSunExt + LWRabs_t.LWRabsWallSun - G1WallSun - HfluxWallSun
    Ycanyon[5] =
        SWRabs_t.SWRabsWallShadeExt + LWRabs_t.LWRabsWallShade - G1WallShade -
        HfluxWallShade

    # Tree energy balance
    if Gemeotry_m.trees > 0
        Ycanyon[6] =
            SWRabs_t.SWRabsTree + LWRabs_t.LWRabsTree - HfluxTree - LEfluxTreeInt -
            LTEfluxTree
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
        Gemoetry_m.hcanyon * (HfluxWallSun + HfluxWallShade) +
        Ctree * 4 * Gemoetry_m.radius_tree * HfluxTree - HfluxCanyon - dS_H_air +
        WasteHeat.SensibleFromVent_Can +
        WasteHeat.SensibleFromAC_Can +
        WasteHeat.SensibleFromHeat_Can

    Ycanyon[10] =
        Cimp * FractionsGround.fimp * LEfluxGroundImp +
        Cbare * FractionsGround.fbare * (LEfluxGroundBarePond + LEfluxGroundBareSoil) +
        Cveg *
        FractionsGround.fveg *
        (LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg) +
        Ctree * 4 * Gemoetry_m.radius_tree * (LEfluxTreeInt + LTEfluxTree) - LEfluxCanyon -
        dS_LE_air +
        WasteHeat.LatentFromVent_Can +
        WasteHeat.LatentFromAC_Can +
        WasteHeat.LatentFromHeat_Can

    return Ycanyon, G2WallSun, G2WallShade, SWRabs_t
end
