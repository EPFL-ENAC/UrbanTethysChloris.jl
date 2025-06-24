"""
    conductive_heat_flux_building_floor(
        Tinground::FT,
        TingroundDamptm1::FT,
        Tingroundtm1::FT,
        ParCalculation::NamedTuple,
        ParThermalBuildFloor::ModelComponents.Parameters.ThermalBuilding{FT}
    ) where {FT<:AbstractFloat}

Calculate conductive heat flux through building floor.

# Arguments
- `Tinground`: Ground temperature [K]
- `TingroundDamptm1`: Damped ground temperature at previous timestep [K]
- `Tingroundtm1`: Ground temperature at previous timestep [K]
- `ParCalculation`: Calculation parameters
- `ParThermalBuildFloor`: Thermal parameters for building floor

# Returns
- `G::FT`: Heat flux through building floor [W/m²]
- `Tdp::FT`: Updated damped temperature [K]
"""
function conductive_heat_flux_building_floor(
    Tinground::FT,
    TingroundDamptm1::FT,
    Tingroundtm1::FT,
    ParCalculation::NamedTuple,
    ParThermalBuildFloor::ModelComponents.Parameters.ThermalBuilding{FT},
) where {FT<:AbstractFloat}
    # Extract parameters
    Ts = Tinground
    Tdptm1 = TingroundDamptm1
    Tstm1 = Tingroundtm1
    dts = ParCalculation.dts
    lan_dry_imp = ParThermalBuildFloor.lan_ground_floor
    cv_s_imp = ParThermalBuildFloor.cv_ground_floor

    # Constants
    tau = FT(86400)  # [s] time constant
    CTt = 2 * sqrt(FT(π) / (lan_dry_imp * cv_s_imp * tau))  # [K m²/J] Total Thermal Capacity Soil

    # Calculate heat flux using soil heat function
    G, Tdp = soil_heat(dts, Ts - FT(273.15), Tstm1 - FT(273.15), Tdptm1 - FT(273.15), CTt)
    Tdp = Tdp + FT(273.15)

    return G, Tdp
end
