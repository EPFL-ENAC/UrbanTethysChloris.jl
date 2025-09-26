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
    _, _, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, _ = Radiation.total_shortwave_absorbed(
        Gemeotry_m,
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

    # Calculate longwave radiation
    _, _, LWRabs_t, _ = Radiation.total_longwave_absorbed(
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
    G1GroundImp, _ = ConductiveHeat.conductive_heat_flux_ground_fr(
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
    G1GroundBare, _ = ConductiveHeat.conductive_heat_flux_ground_vb(
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
    G1GroundVeg, _ = ConductiveHeat.conductive_heat_flux_ground_vb(
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

    # Calculate sensible and latent heat fluxes
    HfluxCanyon, LEfluxCanyon, _, _, _, _ = TurbulentHeat.heat_flux_canyon(
        TemperatureC, Gemeotry_m, MeteoData, ParVegTree, fconvPreCalc, fconv
    )

    # Ground and tree heat fluxes
    HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, _, _, _, _, _, _, _, _, _, _, _, _, LEfluxGroundImp, LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg, LEfluxTreeInt, LTEfluxTree = TurbulentHeat.heat_flux_ground(
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
    HfluxWallSun, HfluxWallShade, _, _, _, _, _, _, _, _, _, _, _, _, cp_atm, rho_atm, L_heat = TurbulentHeat.heat_flux_wall(
        TemperatureC, Gemeotry_m, MeteoData, ParVegTree, ParVegGround, FractionsGround
    )

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

    # Energy balance calculations
    if Gemeotry_m.trees == 0
        SWRabs_t = merge(SWRabs_t, (Tree=zero(FT),))
        LWRabs_t = merge(LWRabs_t, (Tree=zero(FT),))
    end

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

    return Ycanyon,
    G2WallSun, G2WallShade, SWRabs_t, SWRabsWallSunTransmitted,
    SWRabsWallShadeTransmitted
end
