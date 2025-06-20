"""
    conductive_heat_flux_ground_imp(
        TemperatureC::Vector{FT},
        TempVec_ittm::NamedTuple,
        ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParCalculation::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate conductive heat flux through impervious ground.

# Arguments
- `TemperatureC`: Canyon temperatures vector
- `TempVec_ittm`: Temperature vectors at previous time step
- `ParThermalGround`: Thermal parameters for ground
- `ParSoilGround`: Soil parameters for ground
- `ParCalculation`: Calculation parameters

# Returns
- `G1::FT`: Heat flux from surface to concrete interior [W/m²]
- `G2::FT`: Heat flux from concrete interior to deep soil [W/m²]
- `dS::FT`: Energy storage in the ground [J/m²]
"""
function conductive_heat_flux_ground_imp(
    TemperatureC::Vector{FT},
    TempVec_ittm::NamedTuple,
    ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParCalculation::NamedTuple,
) where {FT<:AbstractFloat}
    # Extract inputs
    Ts = TemperatureC[1]
    Tint = TemperatureC[10]
    Tint_tm1 = TempVec_ittm.TGroundIntImp
    Tb = FT(293.15)  # Deep soil temperature [K]

    lan_dry1 = ParThermalGround.lan_dry
    lan_dry2 = ParThermalGround.lan_dry
    dz1 = ParSoilGround.dz1
    dz2 = ParSoilGround.dz2
    cv_s1 = ParThermalGround.cv_s
    cv_s2 = ParThermalGround.cv_s
    dts = ParCalculation.dts

    # Computation of heat fluxes
    G1 = lan_dry1 * (Ts - Tint) / dz1
    G2 = lan_dry2 * (Tint - Tb) / dz2
    dS = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

    return G1, G2, dS
end
