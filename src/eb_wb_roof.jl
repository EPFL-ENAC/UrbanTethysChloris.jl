"""
    eb_solver_roof(
        TemperatureR::Vector{FT},
        TemperatureB::Vector{FT},
        TempVec_ittm::NamedTuple,
        MeteoData::NamedTuple,
        Int_ittm::NamedTuple,
        ExWater_ittm::NamedTuple,
        Vwater_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        SoilPotW_ittm::NamedTuple,
        CiCO2Leaf_ittm::NamedTuple,
        Runon_ittm::NamedTuple,
        Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        PropOpticalRoof::ModelComponents.Parameters.OutdoorOpticalProperties{FT},
        ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        HumidityAtm::NamedTuple,
        Anthropogenic::NamedTuple,
        ParCalculation::NamedTuple,
        BEM_on::Bool,
        RESPreCalc::Bool,
        rsRoofPreCalc::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate energy balance for roof surfaces.

# Arguments
- `TemperatureR`: Roof temperatures [K]
- `TemperatureB`: Building temperatures [K]
- `TempVec_ittm`: Temperature vectors at previous time step
- `MeteoData`: Meteorological data
- `Int_ittm`: Interception at previous time step
- `ExWater_ittm`: Extractable water at previous time step
- `Vwater_ittm`: Soil water volume at previous time step
- `Owater_ittm`: Soil water content at previous time step
- `SoilPotW_ittm`: Soil water potential at previous time step
- `CiCO2Leaf_ittm`: Leaf CO2 concentration at previous time step
- `Runon_ittm`: Runon at previous time step
- `Geometry_m`: Urban geometry parameters
- `FractionsRoof`: Roof surface fractions
- `ParSoilRoof`: Soil parameters for roof
- `PropOpticalRoof`: Optical properties for roof
- `ParThermalRoof`: Thermal properties for roof
- `ParVegRoof`: Vegetation parameters for roof
- `HumidityAtm`: Atmospheric humidity parameters
- `Anthropogenic`: Anthropogenic parameters
- `ParCalculation`: Calculation parameters
- `BEM_on`: Building Energy Model switch
- `RESPreCalc`: Use pre-calculated resistances
- `rsRoofPreCalc`: Pre-calculated resistance parameters

# Returns
- `SWRabsRoofImp`: Absorbed shortwave radiation by impervious roof [W/m²]
- `SWRabsRoofVeg`: Absorbed shortwave radiation by vegetated roof [W/m²]
- `SWRabsTotalRoof`: Total absorbed shortwave radiation by roof [W/m²]
- `SWRoutRoofImp`: Outgoing shortwave radiation from impervious roof [W/m²]
- `SWRoutRoofVeg`: Outgoing shortwave radiation from vegetated roof [W/m²]
- `SWRoutTotalRoof`: Total outgoing shortwave radiation from roof [W/m²]
- `SWRinRoofImp`: Incoming shortwave radiation on impervious roof [W/m²]
- `SWRinRoofVeg`: Incoming shortwave radiation on vegetated roof [W/m²]
- `SWRinTotalRoof`: Total incoming shortwave radiation on roof [W/m²]
- `SWREBRoofImp`: Shortwave radiation energy balance for impervious roof [W/m²]
- `SWREBRoofVeg`: Shortwave radiation energy balance for vegetated roof [W/m²]
- `SWREBTotalRoof`: Total shortwave radiation energy balance for roof [W/m²]
- `LWRabsRoofVeg`: Absorbed longwave radiation by vegetated roof [W/m²]
- `LWRabsRoofImp`: Absorbed longwave radiation by impervious roof [W/m²]
- `LWRabsTotalRoof`: Total absorbed longwave radiation by roof [W/m²]
- `LWRoutRoofVeg`: Outgoing longwave radiation from vegetated roof [W/m²]
- `LWRoutRoofImp`: Outgoing longwave radiation from impervious roof [W/m²]
- `LWRoutTotalRoof`: Total outgoing longwave radiation from roof [W/m²]
- `LWRinRoofImp`: Incoming longwave radiation on impervious roof [W/m²]
- `LWRinRoofVeg`: Incoming longwave radiation on vegetated roof [W/m²]
- `LWRinTotalRoof`: Total incoming longwave radiation on roof [W/m²]
- `LWREBRoofImp`: Longwave radiation energy balance for impervious roof [W/m²]
- `LWREBRoofVeg`: Longwave radiation energy balance for vegetated roof [W/m²]
- `LWREBTotalRoof`: Total longwave radiation energy balance for roof [W/m²]
- `HfluxRoofImp`: Sensible heat flux from impervious roof [W/m²]
- `HfluxRoofVeg`: Sensible heat flux from vegetated roof [W/m²]
- `HfluxRoof`: Total sensible heat flux from roof [W/m²]
- `LEfluxRoofImp`: Latent heat flux from impervious roof [W/m²]
- `LEfluxRoofVegInt`: Interception latent heat flux from vegetated roof [W/m²]
- `LEfluxRoofVegPond`: Pond latent heat flux from vegetated roof [W/m²]
- `LEfluxRoofVegSoil`: Soil latent heat flux from vegetated roof [W/m²]
- `LTEfluxRoofVeg`: Transpiration latent heat flux from vegetated roof [W/m²]
- `LEfluxRoofVeg`: Total latent heat flux from vegetated roof [W/m²]
- `LEfluxRoof`: Total latent heat flux from roof [W/m²]
- `G1RoofImp`: Conductive heat flux through impervious roof [W/m²]
- `G2RoofImp`: Conductive heat flux from impervious roof to building [W/m²]
- `dsRoofImp`: Energy storage change in impervious roof [W/m²]
- `G1RoofVeg`: Conductive heat flux through vegetated roof [W/m²]
- `G2RoofVeg`: Conductive heat flux from vegetated roof to building [W/m²]
- `dsRoofVeg`: Energy storage change in vegetated roof [W/m²]
- `G1Roof`: Total conductive heat flux through roof [W/m²]
- `G2Roof`: Total conductive heat flux from roof to building [W/m²]
- `dsRoof`: Total energy storage change in roof [W/m²]
- `raRooftoAtm`: Aerodynamic resistance from roof to atmosphere [s/m]
- `rb_LRoof`: Leaf boundary layer resistance on roof [s/m]
- `rap_LRoof`: Aerodynamic resistance for pond on roof [s/m]
- `r_soilRoof`: Soil resistance on roof [s/m]
- `rs_sunRoof`: Stomatal resistance for sunlit leaves on roof [s/m]
- `rs_shdRoof`: Stomatal resistance for shaded leaves on roof [s/m]
- `EfluxRoofImp`: Water vapor flux from impervious roof [kg/m²/s]
- `EfluxRoofVegInt`: Water vapor flux from intercepted water on vegetated roof [kg/m²/s]
- `EfluxRoofVegPond`: Water vapor flux from ponded water on vegetated roof [kg/m²/s]
- `EfluxRoofVegSoil`: Water vapor flux from soil on vegetated roof [kg/m²/s]
- `TEfluxRoofVeg`: Transpiration from vegetation on roof [kg/m²/s]
- `EfluxRoofVeg`: Total water vapor flux from vegetated roof [kg/m²/s]
- `EfluxRoof`: Total water vapor flux from roof [kg/m²/s]
- `QRoofImp`: Runoff from impervious roof [mm/dth]
- `QRoofVegDrip`: Drip from vegetated roof [mm/dth]
- `QRoofVegPond`: Ponding on vegetated roof [mm/dth]
- `LkRoofImp`: Leakage from impervious roof [mm/h]
- `LkRoofVeg`: Leakage from vegetated roof [mm/h]
- `LkRoof`: Total roof leakage [mm/h]
- `QRoofVegSoil`: Soil runoff from vegetated roof [mm/dth]
- `RunoffRoofTot`: Total roof runoff [mm/dth]
- `RunonRoofTot`: Total roof runon [mm/dth]
- `IntRoofImp`: Interception on impervious roof [mm]
- `IntRoofVegPlant`: Interception on vegetation [mm]
- `IntRoofVegGround`: Interception on ground below vegetation [mm]
- `dInt_dtRoofImp`: Change in impervious roof interception [mm/dth]
- `dInt_dtRoofVegPlant`: Change in vegetation interception [mm/dth]
- `dInt_dtRoofVegGround`: Change in ground interception [mm/dth]
- `IntRooftot`: Total roof interception [mm]
- `dInt_dtRooftot`: Total change in roof interception [mm/dth]
- `dVRoofSoil_dt`: Change in soil water volume [mm/dth]
- `fRoofVeg`: Infiltration rate on vegetated roof [mm/h]
- `VRoofSoil`: Soil water volume [mm]
- `OwRoofSoil`: Soil water content [-]
- `OSwRoofSoil`: Surface soil water content [-]
- `ExWaterRoof_H`: Extractable water from higher soil layers [mm/h]
- `SoilPotWRoof_H`: Soil water potential in higher soil layers [MPa]
- `ExWaterRoof_L`: Extractable water from lower soil layers [mm/h]
- `SoilPotWRoof_L`: Soil water potential in lower soil layers [MPa]
- `CiCO2LeafRoofVegSun`: Sunlit leaf CO2 concentration [µmol/mol]
- `CiCO2LeafRoofVegShd`: Shaded leaf CO2 concentration [µmol/mol]
- `WBRoofVegInVeg`: Water balance for vegetation interception [mm/dth]
- `WBRoofVegInGround`: Water balance for ground interception [mm/dth]
- `WBRoofVegSoil`: Water balance for soil [mm/dth]
- `EBRoofImp`: Energy balance for impervious roof [W/m²]
- `EBRoofVeg`: Energy balance for vegetated roof [W/m²]
- `Yroof`: Energy balance residuals [W/m²]
- `WBRoofImp`: Water balance for impervious roof [mm/dth]
- `WBRoofVeg`: Water balance for vegetated roof [mm/dth]
- `WBRoofTot`: Total roof water balance [mm/dth]
"""
function eb_wb_roof!(
    model::Model{FT},
    TemperatureR::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    ParCalculation::NamedTuple,
    BEM_on::Bool,
    RESPreCalc::Bool,
    rsRoofPreCalc::NamedTuple,
) where {FT<:AbstractFloat}
    results = eb_wb_roof(
        TemperatureR,
        TemperatureB,
        TempVec_ittm,
        model.forcing.meteorological,
        model.variables.waterflux.Interception,
        model.variables.waterflux.ExWater,
        model.variables.waterflux.Vwater,
        model.variables.waterflux.Owater,
        model.variables.waterflux.SoilPotW,
        model.variables.waterflux.CiCO2Leaf,
        model.variables.waterflux.Runon,
        model.parameters.urbangeometry,
        model.parameters.surfacefractions.roof,
        model.parameters.soil.roof,
        model.parameters.optical.roof,
        model.parameters.thermal.roof,
        model.parameters.vegetation.roof,
        model.forcing.meteorological,
        model.forcing.anthropogenic,
        ParCalculation,
        BEM_on,
        RESPreCalc,
        rsRoofPreCalc,
    )

    update!(model.variables, results, eb_wb_roof_dispatcher)

    return results.G2Roof, results.Yroof
end

function eb_wb_roof(
    TemperatureR::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    Int_ittm::ModelComponents.ModelVariables.Interception{FT},
    ExWater_ittm::ModelComponents.ModelVariables.ExWater{FT,MR,MG},
    Vwater_ittm::ModelComponents.ModelVariables.Vwater{FT,MR,MG},
    Owater_ittm::ModelComponents.ModelVariables.Owater{FT,MR,MG},
    SoilPotW_ittm::ModelComponents.ModelVariables.SoilPotW{FT},
    CiCO2Leaf_ittm::ModelComponents.ModelVariables.CiCO2Leaf{FT},
    Runon_ittm::ModelComponents.ModelVariables.Runon{FT},
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties,
    ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    HumidityAtm::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    Anthropogenic::ModelComponents.ForcingInputs.AnthropogenicInputs{FT,0},
    ParCalculation::NamedTuple,
    BEM_on::Bool,
    RESPreCalc::Bool,
    rsRoofPreCalc::NamedTuple,
) where {FT<:AbstractFloat,MR,MG}
    # TemperatureR(1) = Troof_imp
    # TemperatureR(2) = Troof_veg
    # TemperatureR(3) = Troof_interior_imp
    # TemperatureR(4) = Troof_interior_veg

    # Shortwave radiation
    SWRabs_dir_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_dir  # Absorbed direct shortwave radiation by vegetated roof [W/m²]
    SWRabs_diff_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_diff  # Absorbed diffuse shortwave radiation by vegetated roof [W/m²]
    SWRabsRoofVeg = (1 - PropOpticalRoof.aveg) * (MeteoData.SW_dir + MeteoData.SW_diff)  # Absorbed total shortwave radiation by vegetated roof [W/m²]
    SWRabsRoofImp = (1 - PropOpticalRoof.aimp) * (MeteoData.SW_dir + MeteoData.SW_diff)  # Absorbed total shortwave radiation by impervious roof [W/m²]
    SWRabsTotalRoof =
        SWRabsRoofImp * FractionsRoof.fimp + SWRabsRoofVeg * FractionsRoof.fveg

    SWRoutRoofImp = PropOpticalRoof.aimp * (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRoutRoofVeg = PropOpticalRoof.aveg * (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRoutTotalRoof =
        SWRoutRoofImp * FractionsRoof.fimp + SWRoutRoofVeg * FractionsRoof.fveg

    SWRinRoofImp = (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRinRoofVeg = (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRinTotalRoof = (MeteoData.SW_dir + MeteoData.SW_diff)

    SWREBRoofImp = SWRinRoofImp - SWRoutRoofImp - SWRabsRoofImp
    SWREBRoofVeg = SWRinRoofVeg - SWRoutRoofVeg - SWRabsRoofVeg
    SWREBTotalRoof = SWRinTotalRoof - SWRoutTotalRoof - SWRabsTotalRoof

    # Longwave radiation
    bolzm = FT(5.67e-8)  # Stefan-Boltzmann constant [W*m^-2*K^-4]
    LWRabsRoofVeg =
        MeteoData.LWR - (
            PropOpticalRoof.eveg * bolzm * (TemperatureR[2])^4 +
            (1 - PropOpticalRoof.eveg) * MeteoData.LWR
        )  # Total absorbed longwave radiation by vegetated roof [W/m²]
    LWRabsRoofImp =
        MeteoData.LWR - (
            PropOpticalRoof.eimp * bolzm * (TemperatureR[1])^4 +
            (1 - PropOpticalRoof.eimp) * MeteoData.LWR
        )  # Total absorbed longwave radiation by impervious roof [W/m²]
    LWRabsTotalRoof =
        LWRabsRoofImp * FractionsRoof.fimp + LWRabsRoofVeg * FractionsRoof.fveg

    LWRoutRoofVeg =
        PropOpticalRoof.eveg * bolzm * (TemperatureR[2])^4 +
        (1 - PropOpticalRoof.eveg) * MeteoData.LWR
    LWRoutRoofImp =
        PropOpticalRoof.eimp * bolzm * (TemperatureR[1])^4 +
        (1 - PropOpticalRoof.eimp) * MeteoData.LWR
    LWRoutTotalRoof =
        LWRoutRoofImp * FractionsRoof.fimp + LWRoutRoofVeg * FractionsRoof.fveg

    LWRinRoofImp = MeteoData.LWR
    LWRinRoofVeg = MeteoData.LWR
    LWRinTotalRoof = MeteoData.LWR

    LWREBRoofImp = LWRinRoofImp - LWRoutRoofImp - LWRabsRoofImp
    LWREBRoofVeg = LWRinRoofVeg - LWRoutRoofVeg - LWRabsRoofVeg
    LWREBTotalRoof = LWRinTotalRoof - LWRoutTotalRoof - LWRabsTotalRoof

    # Sensible and latent heat
    HfluxRoofImp, HfluxRoofVeg, EfluxRoofImp, EfluxRoofVegInt, EfluxRoofVegPond, EfluxRoofVegSoil, TEfluxRoofVeg, LEfluxRoofImp, LEfluxRoofVegInt, LEfluxRoofVegPond, LEfluxRoofVegSoil, LTEfluxRoofVeg, CiCO2LeafRoofVegSun, CiCO2LeafRoofVegShd, raRooftoAtm, rb_LRoof, rap_LRoof, r_soilRoof, rs_sunRoof, rs_shdRoof = TurbulentHeat.heat_flux_roof(
        TemperatureR,
        TempVec_ittm,
        MeteoData,
        HumidityAtm,
        ParVegRoof,
        FractionsRoof,
        Geometry_m,
        ParSoilRoof,
        ParCalculation,
        SoilPotW_ittm,
        Owater_ittm,
        Vwater_ittm,
        ExWater_ittm,
        Int_ittm,
        CiCO2Leaf_ittm,
        SWRabs_dir_veg,
        SWRabs_diff_veg,
        RESPreCalc,
        rsRoofPreCalc,
    )

    HfluxRoof = HfluxRoofImp * FractionsRoof.fimp + HfluxRoofVeg * FractionsRoof.fveg
    LEfluxRoofVeg =
        LEfluxRoofVegInt + LEfluxRoofVegPond + LEfluxRoofVegSoil + LTEfluxRoofVeg
    LEfluxRoof = LEfluxRoofImp * FractionsRoof.fimp + LEfluxRoofVeg * FractionsRoof.fveg
    EfluxRoofVeg = EfluxRoofVegInt + EfluxRoofVegPond + EfluxRoofVegSoil + TEfluxRoofVeg
    EfluxRoof = EfluxRoofImp * FractionsRoof.fimp + EfluxRoofVeg * FractionsRoof.fveg

    # Conductive heat fluxes roof
    # Impervious conductive heat flux
    G1RoofImp, G2RoofImp, dsRoofImp = ConductiveHeat.conductive_heat_flux_roof_imp(
        TemperatureR,
        TemperatureB,
        TempVec_ittm,
        Anthropogenic,
        ParThermalRoof,
        ParSoilRoof,
        ParCalculation,
        BEM_on,
    )

    # Conductive heat flux of green roof
    G1RoofVeg, G2RoofVeg, dsRoofVeg = ConductiveHeat.conductive_heat_flux_green_roof(
        TemperatureR,
        TemperatureB,
        TempVec_ittm,
        Anthropogenic,
        Owater_ittm,
        ParVegRoof,
        ParSoilRoof,
        ParThermalRoof,
        ParCalculation,
        BEM_on,
    )

    G1Roof = G1RoofImp * FractionsRoof.fimp + G1RoofVeg * FractionsRoof.fveg
    G2Roof = G2RoofImp * FractionsRoof.fimp + G2RoofVeg * FractionsRoof.fveg
    dsRoof = dsRoofImp * FractionsRoof.fimp + dsRoofVeg * FractionsRoof.fveg

    # Energy balance
    Yroof = zeros(FT, 4)

    if FractionsRoof.fimp > 0
        Yroof[1] = SWRabsRoofImp + LWRabsRoofImp - HfluxRoofImp - G1RoofImp - LEfluxRoofImp  # Energy budget impervious roof
        Yroof[3] = G1RoofImp - G2RoofImp - dsRoofImp  # Energy budget concrete mass roof
    else
        Yroof[1] = TemperatureR[1] - FT(273.15)
        Yroof[3] = TemperatureR[3] - FT(273.15)
    end

    if FractionsRoof.fveg > 0
        Yroof[2] =
            SWRabsRoofVeg + LWRabsRoofVeg - HfluxRoofVeg - G1RoofVeg - LEfluxRoofVegInt -
            LEfluxRoofVegPond - LEfluxRoofVegSoil - LTEfluxRoofVeg  # Energy budget vegetated roof
        Yroof[4] = G1RoofVeg - G2RoofVeg - dsRoofVeg  # Energy budget concrete mass roof
    else
        Yroof[2] = TemperatureR[2] - FT(273.15)
        Yroof[4] = TemperatureR[4] - FT(273.15)
    end

    EBRoofImp =
        SWRabsRoofImp + LWRabsRoofImp - HfluxRoofImp - G2RoofImp - dsRoofImp - LEfluxRoofImp
    EBRoofVeg =
        SWRabsRoofVeg + LWRabsRoofVeg - HfluxRoofVeg - G2RoofVeg - dsRoofVeg -
        LEfluxRoofVegInt - LEfluxRoofVegPond - LEfluxRoofVegSoil - LTEfluxRoofVeg

    # Water fluxes and water balance
    QRoofImp, IntRoofImp, dInt_dtRoofImp, LkRoofImp, QRoofVegDrip, IntRoofVegPlant, dInt_dtRoofVegPlant, QRoofVegPond, IntRoofVegGround, dInt_dtRoofVegGround, dVRoofSoil_dt, fRoofVeg, VRoofSoil, OwRoofSoil, OSwRoofSoil, LkRoofVeg, SoilPotWRoof_L, ExWaterRoof_L, QRoofVegSoil, TEfluxRoofVeg, EfluxRoofVegSoil, RunoffRoofTot, RunonRoofTot, WBRoofVegInVeg, WBRoofVegInGround, WBRoofVegSoil, WBRoofImp, WBRoofVeg, WBRoofTot = Water.water_roof(
        EfluxRoofImp,
        EfluxRoofVegInt,
        EfluxRoofVegPond,
        EfluxRoofVegSoil,
        TEfluxRoofVeg,
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        FractionsRoof,
        ParSoilRoof,
        ParCalculation,
        ParVegRoof,
        Anthropogenic,
    )

    LkRoof = LkRoofImp * FractionsRoof.fimp + LkRoofVeg * FractionsRoof.fveg
    IntRooftot =
        IntRoofImp * FractionsRoof.fimp +
        (IntRoofVegPlant + IntRoofVegGround) * FractionsRoof.fveg
    dInt_dtRooftot =
        dInt_dtRoofImp * FractionsRoof.fimp +
        (dInt_dtRoofVegPlant + dInt_dtRoofVegGround) * FractionsRoof.fveg
    ExWaterRoof_H = fill(FT(NaN), ParSoilRoof.ms)
    SoilPotWRoof_H = FT(NaN)

    return (;
        SWRabsRoofImp,
        SWRabsRoofVeg,
        SWRabsTotalRoof,
        SWRoutRoofImp,
        SWRoutRoofVeg,
        SWRoutTotalRoof,
        SWRinRoofImp,
        SWRinRoofVeg,
        SWRinTotalRoof,
        SWREBRoofImp,
        SWREBRoofVeg,
        SWREBTotalRoof,
        LWRabsRoofVeg,
        LWRabsRoofImp,
        LWRabsTotalRoof,
        LWRoutRoofVeg,
        LWRoutRoofImp,
        LWRoutTotalRoof,
        LWRinRoofImp,
        LWRinRoofVeg,
        LWRinTotalRoof,
        LWREBRoofImp,
        LWREBRoofVeg,
        LWREBTotalRoof,
        HfluxRoofImp,
        HfluxRoofVeg,
        HfluxRoof,
        LEfluxRoofImp,
        LEfluxRoofVegInt,
        LEfluxRoofVegPond,
        LEfluxRoofVegSoil,
        LTEfluxRoofVeg,
        LEfluxRoofVeg,
        LEfluxRoof,
        G1RoofImp,
        G2RoofImp,
        dsRoofImp,
        G1RoofVeg,
        G2RoofVeg,
        dsRoofVeg,
        G1Roof,
        G2Roof,
        dsRoof,
        raRooftoAtm,
        rb_LRoof,
        rap_LRoof,
        r_soilRoof,
        rs_sunRoof,
        rs_shdRoof,
        EfluxRoofImp,
        EfluxRoofVegInt,
        EfluxRoofVegPond,
        EfluxRoofVegSoil,
        TEfluxRoofVeg,
        EfluxRoofVeg,
        EfluxRoof,
        QRoofImp,
        QRoofVegDrip,
        QRoofVegPond,
        LkRoofImp,
        LkRoofVeg,
        LkRoof,
        QRoofVegSoil,
        RunoffRoofTot,
        RunonRoofTot,
        IntRoofImp,
        IntRoofVegPlant,
        IntRoofVegGround,
        dInt_dtRoofImp,
        dInt_dtRoofVegPlant,
        dInt_dtRoofVegGround,
        IntRooftot,
        dInt_dtRooftot,
        dVRoofSoilVeg_dt=dVRoofSoil_dt,
        fRoofVeg,
        VRoofSoilVeg=VRoofSoil,
        OwRoofSoilVeg=OwRoofSoil,
        OSwRoofSoilVeg=OSwRoofSoil,
        ExWaterRoofVeg_H=ExWaterRoof_H,
        SoilPotWRoofVeg_H=SoilPotWRoof_H,
        SoilPotWRoofVeg_L=SoilPotWRoof_L,
        ExWaterRoofVeg_L=ExWaterRoof_L,
        CiCO2LeafRoofVegSun,
        CiCO2LeafRoofVegShd,
        WBRoofVegInVeg,
        WBRoofVegInGround,
        WBRoofVegSoil,
        EBRoofImp,
        EBRoofVeg,
        Yroof,
        WBRoofImp,
        WBRoofVeg,
        WBRoofTot,
    )
end
function eb_wb_roof(
    TemperatureR::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::NamedTuple,
    MeteoData::NamedTuple,
    Int_ittm::NamedTuple,
    ExWater_ittm::NamedTuple,
    Vwater_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    SoilPotW_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    Runon_ittm::NamedTuple,
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties,
    ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    HumidityAtm::NamedTuple,
    Anthropogenic::NamedTuple,
    ParCalculation::NamedTuple,
    BEM_on::Bool,
    RESPreCalc::Bool,
    rsRoofPreCalc::NamedTuple,
) where {FT<:AbstractFloat}
    # TemperatureR(1) = Troof_imp
    # TemperatureR(2) = Troof_veg
    # TemperatureR(3) = Troof_interior_imp
    # TemperatureR(4) = Troof_interior_veg

    # Shortwave radiation
    SWRabs_dir_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_dir  # Absorbed direct shortwave radiation by vegetated roof [W/m²]
    SWRabs_diff_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_diff  # Absorbed diffuse shortwave radiation by vegetated roof [W/m²]
    SWRabsRoofVeg = (1 - PropOpticalRoof.aveg) * (MeteoData.SW_dir + MeteoData.SW_diff)  # Absorbed total shortwave radiation by vegetated roof [W/m²]
    SWRabsRoofImp = (1 - PropOpticalRoof.aimp) * (MeteoData.SW_dir + MeteoData.SW_diff)  # Absorbed total shortwave radiation by impervious roof [W/m²]
    SWRabsTotalRoof =
        SWRabsRoofImp * FractionsRoof.fimp + SWRabsRoofVeg * FractionsRoof.fveg

    SWRoutRoofImp = PropOpticalRoof.aimp * (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRoutRoofVeg = PropOpticalRoof.aveg * (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRoutTotalRoof =
        SWRoutRoofImp * FractionsRoof.fimp + SWRoutRoofVeg * FractionsRoof.fveg

    SWRinRoofImp = (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRinRoofVeg = (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRinTotalRoof = (MeteoData.SW_dir + MeteoData.SW_diff)

    SWREBRoofImp = SWRinRoofImp - SWRoutRoofImp - SWRabsRoofImp
    SWREBRoofVeg = SWRinRoofVeg - SWRoutRoofVeg - SWRabsRoofVeg
    SWREBTotalRoof = SWRinTotalRoof - SWRoutTotalRoof - SWRabsTotalRoof

    # Longwave radiation
    bolzm = FT(5.67e-8)  # Stefan-Boltzmann constant [W*m^-2*K^-4]
    LWRabsRoofVeg =
        MeteoData.LWR - (
            PropOpticalRoof.eveg * bolzm * (TemperatureR[2])^4 +
            (1 - PropOpticalRoof.eveg) * MeteoData.LWR
        )  # Total absorbed longwave radiation by vegetated roof [W/m²]
    LWRabsRoofImp =
        MeteoData.LWR - (
            PropOpticalRoof.eimp * bolzm * (TemperatureR[1])^4 +
            (1 - PropOpticalRoof.eimp) * MeteoData.LWR
        )  # Total absorbed longwave radiation by impervious roof [W/m²]
    LWRabsTotalRoof =
        LWRabsRoofImp * FractionsRoof.fimp + LWRabsRoofVeg * FractionsRoof.fveg

    LWRoutRoofVeg =
        PropOpticalRoof.eveg * bolzm * (TemperatureR[2])^4 +
        (1 - PropOpticalRoof.eveg) * MeteoData.LWR
    LWRoutRoofImp =
        PropOpticalRoof.eimp * bolzm * (TemperatureR[1])^4 +
        (1 - PropOpticalRoof.eimp) * MeteoData.LWR
    LWRoutTotalRoof =
        LWRoutRoofImp * FractionsRoof.fimp + LWRoutRoofVeg * FractionsRoof.fveg

    LWRinRoofImp = MeteoData.LWR
    LWRinRoofVeg = MeteoData.LWR
    LWRinTotalRoof = MeteoData.LWR

    LWREBRoofImp = LWRinRoofImp - LWRoutRoofImp - LWRabsRoofImp
    LWREBRoofVeg = LWRinRoofVeg - LWRoutRoofVeg - LWRabsRoofVeg
    LWREBTotalRoof = LWRinTotalRoof - LWRoutTotalRoof - LWRabsTotalRoof

    # Sensible and latent heat
    HfluxRoofImp, HfluxRoofVeg, EfluxRoofImp, EfluxRoofVegInt, EfluxRoofVegPond, EfluxRoofVegSoil, TEfluxRoofVeg, LEfluxRoofImp, LEfluxRoofVegInt, LEfluxRoofVegPond, LEfluxRoofVegSoil, LTEfluxRoofVeg, CiCO2LeafRoofVegSun, CiCO2LeafRoofVegShd, raRooftoAtm, rb_LRoof, rap_LRoof, r_soilRoof, rs_sunRoof, rs_shdRoof = TurbulentHeat.heat_flux_roof(
        TemperatureR,
        TempVec_ittm,
        MeteoData,
        HumidityAtm,
        ParVegRoof,
        FractionsRoof,
        Geometry_m,
        ParSoilRoof,
        ParCalculation,
        SoilPotW_ittm,
        Owater_ittm,
        Vwater_ittm,
        ExWater_ittm,
        Int_ittm,
        CiCO2Leaf_ittm,
        SWRabs_dir_veg,
        SWRabs_diff_veg,
        RESPreCalc,
        rsRoofPreCalc,
    )

    HfluxRoof = HfluxRoofImp * FractionsRoof.fimp + HfluxRoofVeg * FractionsRoof.fveg
    LEfluxRoofVeg =
        LEfluxRoofVegInt + LEfluxRoofVegPond + LEfluxRoofVegSoil + LTEfluxRoofVeg
    LEfluxRoof = LEfluxRoofImp * FractionsRoof.fimp + LEfluxRoofVeg * FractionsRoof.fveg
    EfluxRoofVeg = EfluxRoofVegInt + EfluxRoofVegPond + EfluxRoofVegSoil + TEfluxRoofVeg
    EfluxRoof = EfluxRoofImp * FractionsRoof.fimp + EfluxRoofVeg * FractionsRoof.fveg

    # Conductive heat fluxes roof
    # Impervious conductive heat flux
    G1RoofImp, G2RoofImp, dsRoofImp = ConductiveHeat.conductive_heat_flux_roof_imp(
        TemperatureR,
        TemperatureB,
        TempVec_ittm,
        Anthropogenic,
        ParThermalRoof,
        ParSoilRoof,
        ParCalculation,
        BEM_on,
    )

    # Conductive heat flux of green roof
    G1RoofVeg, G2RoofVeg, dsRoofVeg = ConductiveHeat.conductive_heat_flux_green_roof(
        TemperatureR,
        TemperatureB,
        TempVec_ittm,
        Anthropogenic,
        Owater_ittm,
        ParVegRoof,
        ParSoilRoof,
        ParThermalRoof,
        ParCalculation,
        BEM_on,
    )

    G1Roof = G1RoofImp * FractionsRoof.fimp + G1RoofVeg * FractionsRoof.fveg
    G2Roof = G2RoofImp * FractionsRoof.fimp + G2RoofVeg * FractionsRoof.fveg
    dsRoof = dsRoofImp * FractionsRoof.fimp + dsRoofVeg * FractionsRoof.fveg

    # Energy balance
    Yroof = zeros(FT, 4)

    if FractionsRoof.fimp > 0
        Yroof[1] = SWRabsRoofImp + LWRabsRoofImp - HfluxRoofImp - G1RoofImp - LEfluxRoofImp  # Energy budget impervious roof
        Yroof[3] = G1RoofImp - G2RoofImp - dsRoofImp  # Energy budget concrete mass roof
    else
        Yroof[1] = TemperatureR[1] - FT(273.15)
        Yroof[3] = TemperatureR[3] - FT(273.15)
    end

    if FractionsRoof.fveg > 0
        Yroof[2] =
            SWRabsRoofVeg + LWRabsRoofVeg - HfluxRoofVeg - G1RoofVeg - LEfluxRoofVegInt -
            LEfluxRoofVegPond - LEfluxRoofVegSoil - LTEfluxRoofVeg  # Energy budget vegetated roof
        Yroof[4] = G1RoofVeg - G2RoofVeg - dsRoofVeg  # Energy budget concrete mass roof
    else
        Yroof[2] = TemperatureR[2] - FT(273.15)
        Yroof[4] = TemperatureR[4] - FT(273.15)
    end

    EBRoofImp =
        SWRabsRoofImp + LWRabsRoofImp - HfluxRoofImp - G2RoofImp - dsRoofImp - LEfluxRoofImp
    EBRoofVeg =
        SWRabsRoofVeg + LWRabsRoofVeg - HfluxRoofVeg - G2RoofVeg - dsRoofVeg -
        LEfluxRoofVegInt - LEfluxRoofVegPond - LEfluxRoofVegSoil - LTEfluxRoofVeg

    # Water fluxes and water balance
    QRoofImp, IntRoofImp, dInt_dtRoofImp, LkRoofImp, QRoofVegDrip, IntRoofVegPlant, dInt_dtRoofVegPlant, QRoofVegPond, IntRoofVegGround, dInt_dtRoofVegGround, dVRoofSoil_dt, fRoofVeg, VRoofSoil, OwRoofSoil, OSwRoofSoil, LkRoofVeg, SoilPotWRoof_L, ExWaterRoof_L, QRoofVegSoil, TEfluxRoofVeg, EfluxRoofVegSoil, RunoffRoofTot, RunonRoofTot, WBRoofVegInVeg, WBRoofVegInGround, WBRoofVegSoil, WBRoofImp, WBRoofVeg, WBRoofTot = Water.water_roof(
        EfluxRoofImp,
        EfluxRoofVegInt,
        EfluxRoofVegPond,
        EfluxRoofVegSoil,
        TEfluxRoofVeg,
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        FractionsRoof,
        ParSoilRoof,
        ParCalculation,
        ParVegRoof,
        Anthropogenic,
    )

    LkRoof = LkRoofImp * FractionsRoof.fimp + LkRoofVeg * FractionsRoof.fveg
    IntRooftot =
        IntRoofImp * FractionsRoof.fimp +
        (IntRoofVegPlant + IntRoofVegGround) * FractionsRoof.fveg
    dInt_dtRooftot =
        dInt_dtRoofImp * FractionsRoof.fimp +
        (dInt_dtRoofVegPlant + dInt_dtRoofVegGround) * FractionsRoof.fveg
    ExWaterRoof_H = fill(FT(NaN), ParSoilRoof.ms)
    SoilPotWRoof_H = fill(FT(NaN), 1)

    return (;
        SWRabsRoofImp,
        SWRabsRoofVeg,
        SWRabsTotalRoof,
        SWRoutRoofImp,
        SWRoutRoofVeg,
        SWRoutTotalRoof,
        SWRinRoofImp,
        SWRinRoofVeg,
        SWRinTotalRoof,
        SWREBRoofImp,
        SWREBRoofVeg,
        SWREBTotalRoof,
        LWRabsRoofVeg,
        LWRabsRoofImp,
        LWRabsTotalRoof,
        LWRoutRoofVeg,
        LWRoutRoofImp,
        LWRoutTotalRoof,
        LWRinRoofImp,
        LWRinRoofVeg,
        LWRinTotalRoof,
        LWREBRoofImp,
        LWREBRoofVeg,
        LWREBTotalRoof,
        HfluxRoofImp,
        HfluxRoofVeg,
        HfluxRoof,
        LEfluxRoofImp,
        LEfluxRoofVegInt,
        LEfluxRoofVegPond,
        LEfluxRoofVegSoil,
        LTEfluxRoofVeg,
        LEfluxRoofVeg,
        LEfluxRoof,
        G1RoofImp,
        G2RoofImp,
        dsRoofImp,
        G1RoofVeg,
        G2RoofVeg,
        dsRoofVeg,
        G1Roof,
        G2Roof,
        dsRoof,
        raRooftoAtm,
        rb_LRoof,
        rap_LRoof,
        r_soilRoof,
        rs_sunRoof,
        rs_shdRoof,
        EfluxRoofImp,
        EfluxRoofVegInt,
        EfluxRoofVegPond,
        EfluxRoofVegSoil,
        TEfluxRoofVeg,
        EfluxRoofVeg,
        EfluxRoof,
        QRoofImp,
        QRoofVegDrip,
        QRoofVegPond,
        LkRoofImp,
        LkRoofVeg,
        LkRoof,
        QRoofVegSoil,
        RunoffRoofTot,
        RunonRoofTot,
        IntRoofImp,
        IntRoofVegPlant,
        IntRoofVegGround,
        dInt_dtRoofImp,
        dInt_dtRoofVegPlant,
        dInt_dtRoofVegGround,
        IntRooftot,
        dInt_dtRooftot,
        dVRoofSoil_dt,
        fRoofVeg,
        VRoofSoil,
        OwRoofSoil,
        OSwRoofSoil,
        ExWaterRoof_H,
        SoilPotWRoof_H,
        SoilPotWRoof_L,
        ExWaterRoof_L,
        CiCO2LeafRoofVegSun,
        CiCO2LeafRoofVegShd,
        WBRoofVegInVeg,
        WBRoofVegInGround,
        WBRoofVegSoil,
        EBRoofImp,
        EBRoofVeg,
        Yroof,
        WBRoofImp,
        WBRoofVeg,
        WBRoofTot,
    )
end
