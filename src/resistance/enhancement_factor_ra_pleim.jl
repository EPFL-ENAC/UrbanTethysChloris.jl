"""
    enhancement_factor_ra_pleim(
        ra::FT,
        zom::FT,
        zoh::FT,
        disp_h::FT,
        zatm::FT,
        Ws::FT,
        hPBL::FT
    ) where {FT<:AbstractFloat}

Calculate enhancement factor for aerodynamic resistance according to Pleim et al. 2007.

# Arguments
- `ra`: Aerodynamic resistance [s/m]
- `zom`: Momentum roughness length [m]
- `zoh`: Heat roughness length [m]
- `disp_h`: Displacement height [m]
- `zatm`: Atmospheric measurement height [m]
- `Ws`: Wind speed [m/s]
- `hPBL`: Planetary boundary layer height [m]

# Returns
- `fconv`: Fraction of convective transport [-]
- `ra_enhanced`: Enhanced aerodynamic resistance [s/m]
- `ra`: Original aerodynamic resistance [s/m]
- `LAN`: Monin-Obhukov length [m]
"""
function enhancement_factor_ra_pleim(
    ra::FT, zom::FT, zoh::FT, disp_h::FT, zatm::FT, Ws::FT, hPBL::FT
) where {FT<:AbstractFloat}
    # Backcalculate Monin-Obhukov Length
    LAN, _ = backcalculate_obhukov_length(ra, zom, zoh, disp_h, zatm, Ws)

    # Constants from Holtslag et al. 1993
    a = FT(7.2)
    k = FT(0.4)  # von Karman constant

    # Fraction of convective (non-local transport) according to Pleim et al. 2007
    fconv = (1 + k^(-2/3) / (FT(0.1) * a) * (-hPBL/LAN)^(-1/3))^(-1)
    ra_enhanced = ra * (1 - fconv)

    return fconv, ra_enhanced, ra, LAN
end
