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
        Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
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
- `Int_ittm`: Previous timestep interception values
- `ExWater_ittm`: Previous timestep extractable water values
- `Vwater_ittm`: Previous timestep soil water volume values
- `Owater_ittm`: Previous timestep soil water content values
- `SoilPotW_ittm`: Previous timestep soil water potential values
- `CiCO2Leaf_ittm`: Previous timestep leaf CO2 concentration values
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
- `Yroof`: Energy balance residuals [W/m²]
- `G2Roof`: Conductive heat flux through roof [W/m²]
"""
function eb_solver_roof(
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
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    HumidityAtm::NamedTuple,
    Anthropogenic::NamedTuple,
    ParCalculation::NamedTuple,
    BEM_on::Bool,
    RESPreCalc::Bool,
    rsRoofPreCalc::NamedTuple,
) where {FT<:AbstractFloat}
    # Shortwave radiation
    SWRabs_dir_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_dir
    SWRabs_diff_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_diff
    SWRabs_roofveg = (1 - PropOpticalRoof.aveg) * (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRabs_roofimp = (1 - PropOpticalRoof.aimp) * (MeteoData.SW_dir + MeteoData.SW_diff)

    SWRout_roofveg = PropOpticalRoof.aveg * (MeteoData.SW_dir + MeteoData.SW_diff)
    SWRout_roofimp = PropOpticalRoof.aimp * (MeteoData.SW_dir + MeteoData.SW_diff)

    # Longwave radiation
    bolzm = FT(5.67e-8)  # Stefan-Boltzmann constant [W*m^-2*K^-4]
    LWRabs_roofveg =
        MeteoData.LWR - (
            PropOpticalRoof.eveg * bolzm * (TemperatureR[2])^4 +
            (1 - PropOpticalRoof.eveg) * MeteoData.LWR
        )
    LWRabs_roofimp =
        MeteoData.LWR - (
            PropOpticalRoof.eimp * bolzm * (TemperatureR[1])^4 +
            (1 - PropOpticalRoof.eimp) * MeteoData.LWR
        )

    LWRout_roofveg =
        PropOpticalRoof.eveg * bolzm * (TemperatureR[2])^4 +
        (1 - PropOpticalRoof.eveg) * MeteoData.LWR
    LWRout_roofimp =
        PropOpticalRoof.eimp * bolzm * (TemperatureR[1])^4 +
        (1 - PropOpticalRoof.eimp) * MeteoData.LWR

    # Sensible and latent heat
    Hroof_imp, Hroof_veg, Eroof_imp, Eroof_veg, Eroof_ground, Eroof_soil, TEroof_veg, LEroof_imp, LEroof_veg, LEroof_ground, LEroof_soil, LTEroof_veg, Ci_sun_roof, Ci_shd_roof, ra, rb_L, rap_L, r_soil, rs_sun, rs_shd = TurbulentHeat.heat_flux_roof(
        TemperatureR,
        TempVec_ittm,
        MeteoData,
        HumidityAtm,
        ParVegRoof,
        FractionsRoof,
        Gemeotry_m,
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

    # Conductive heat fluxes
    # Impervious conductive heat flux
    G1_roofimp, G2_roofimp, dS_roofimp = ConductiveHeat.conductive_heat_flux_roof_imp(
        TemperatureR,
        TemperatureB,
        TempVec_ittm,
        Anthropogenic,
        ParThermalRoof,
        ParSoilRoof,
        ParCalculation,
        BEM_on,
    )

    # Green roof conductive heat flux
    G1_roofveg, G2_roofveg, dS_roofveg = ConductiveHeat.conductive_heat_flux_green_roof(
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

    G2Roof = FractionsRoof.fimp * G2_roofimp + FractionsRoof.fveg * G2_roofveg

    # Energy balance residuals
    Yroof = zeros(FT, 4)

    if FractionsRoof.fimp > 0
        Yroof[1] = SWRabs_roofimp + LWRabs_roofimp - Hroof_imp - G1_roofimp - LEroof_imp
        Yroof[3] = G1_roofimp - G2_roofimp - dS_roofimp
    else
        Yroof[1] = TemperatureR[1] - FT(273.15)
        Yroof[3] = TemperatureR[3] - FT(273.15)
    end

    if FractionsRoof.fveg > 0
        Yroof[2] =
            SWRabs_roofveg + LWRabs_roofveg - Hroof_veg - G1_roofveg - LEroof_veg -
            LEroof_ground - LEroof_soil - LTEroof_veg
        Yroof[4] = G1_roofveg - G2_roofveg - dS_roofveg
    else
        Yroof[2] = TemperatureR[2] - FT(273.15)
        Yroof[4] = TemperatureR[4] - FT(273.15)
    end

    return Yroof, G2Roof
end
