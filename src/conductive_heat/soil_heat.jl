"""
    soil_heat(dt::Int, Ts::FT, Tstm1::FT, Tdptm1::FT, CTt::FT)


Compute soil heat flux and depth temperature using force-restore method.

# Arguments
- `dt`: time step [s]
- `Ts`: surface temperature [°C]
- `Tstm1`: surface temperature at previous step [°C]
- `Tdptm1`: depth temperature at previous step [°C]
- `CTt`: total thermal capacity of soil or water [K m²/J]

# Returns
- `G::FT`: soil heat flux [W/m²]
- `Tdp::FT`: depth temperature [°C]
"""
function soil_heat(
    dt::Int, Ts::FT, Tstm1::FT, Tdptm1::FT, CTt::FT
) where {FT<:AbstractFloat}

    # Time constant [s]
    τ = FT(86400.0)

    # Temperature variation [°C]
    ΔTs = Ts - Tstm1

    # # Depth temperature calculation using Force Restore (Noilhan & Planton, 1989)
    Tdp = (1.0 / (1.0 + dt / τ)) * (Tdptm1 + (dt / τ) * Ts)

    # Soil heat flux calculation [W/m²]
    G = (1.0 / CTt) * (2π * (Ts - Tdp) / τ + ΔTs / dt)

    return G, Tdp
end
