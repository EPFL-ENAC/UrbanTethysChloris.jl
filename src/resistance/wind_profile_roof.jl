"""
    wind_profile_roof(
        Hcan::FT,
        Zatm::FT,
        uatm::FT,
        zom::FT,
        disp_h::FT,
        hveg::FT,
        Zp::FT
    ) where {FT<:AbstractFloat}

Calculate wind profile above roof surface.

# Arguments
- `Hcan`: Canyon height [m]
- `Zatm`: Atmospheric reference height [m]
- `uatm`: Wind speed at atmospheric reference height [m/s]
- `zom`: Roughness length of roof cover [m]
- `disp_h`: Displacement height [m]
- `hveg`: Vegetation height [m]
- `Zp`: Height of interest [m]

# Returns
- `u_Zp`: Wind speed at height Zp [m/s]
- `u_Hveg`: Wind speed at vegetation top [m/s]
"""
function wind_profile_roof(
    Hcan::FT, Zatm::FT, uatm::FT, zom::FT, disp_h::FT, hveg::FT, Zp::FT
) where {FT<:AbstractFloat}
    k = FT(0.4)  # von Karman constant
    Zatm = Zatm - Hcan
    us_atm = k * uatm / log(Zatm/zom)  # Friction velocity [m/s]
    u_Zp = (us_atm/k) * log((Zp - disp_h)/zom)  # Wind speed at height Zp [m/s]
    u_Hveg = (us_atm/k) * log((hveg - disp_h)/zom)  # Wind speed at vegetation top [m/s]

    return u_Zp, u_Hveg
end
