"""
    backcalculate_obhukov_length(
        ra::FT,
        zom::FT,
        zoh::FT,
        disp_h::FT,
        zatm::FT,
        u::FT
    ) where {FT<:AbstractFloat}

Calculate Monin-Obhukov length based on previously calculated ra to compute enhancement term due to non-local transport in the urban canopy.

# Arguments
- `ra`: Aerodynamic resistance [s/m]
- `zom`: Momentum roughness length [m]
- `zoh`: Heat roughness length [m]
- `disp_h`: Displacement height [m]
- `zatm`: Atmospheric measurement height [m]
- `u`: Wind speed [m/s]

# Returns
- `LAN`: Monin-Obhukov length [m]
- `Diff_ra`: Difference in aerodynamic resistance [s/m]
"""
function backcalculate_obhukov_length(
    ra::FT, zom::FT, zoh::FT, disp_h::FT, zatm::FT, u::FT
) where {FT<:AbstractFloat}
    z = zatm - disp_h
    LANp = FT(-2)  # Initial guess for unstable conditions

    # Find root of solve_obhukov_length
    LAN = find_zero(
        LAN -> solve_obhukov_length(LAN, ra, zom, zoh, z, u), LANp; abstol=FT(1e-6)
    )
    Diff_ra = solve_obhukov_length(LAN, ra, zom, zoh, z, u)

    return LAN, Diff_ra
end

"""
    businger_stability_functions(y::FT) where {FT<:AbstractFloat}

Calculate Businger stability functions for momentum and heat transfer.

# Arguments
- `y`: Stability parameter [-]

# Returns
- `Fih`: Heat stability function [-]
- `Fim`: Momentum stability function [-]
"""
function businger_stability_functions(y::FT) where {FT<:AbstractFloat}
    if y > 0
        # Stable condition
        a = FT(1)
        b = FT(0.667)
        c = FT(5)
        d = FT(0.35)

        Fim = -(a * y + b * (y - c/d) * exp(-d * y) + b * c/d)
        Fih = -(((1 + 2 * a * y/3)^1.5) + b * (y - c/d) * exp(-d * y) + (b * c/d - 1))
    else
        # Unstable condition
        G = FT(16)  # Dyer (1974)
        x = (1 - G * y)^(FT(0.25))
        Fim = log((FT(0.5) * (1 + x^2)) * ((FT(0.5) * (1 + x))^2)) - 2 * atan(x) + π/2
        Fih = 2 * log(FT(0.5) * (1 + x^2))
    end
    return Fih, Fim
end

"""
    solve_obhukov_length(LAN::FT, ra::FT, zom::FT, zoh::FT, z::FT, u::FT) where {FT<:AbstractFloat}

Calculate residual for Monin-Obhukov length solver.

# Arguments
- `LAN`: Monin-Obhukov length [m]
- `ra`: Aerodynamic resistance [s/m]
- `zom`: Momentum roughness length [m]
- `zoh`: Heat roughness length [m]
- `z`: Reference height [m]
- `u`: Wind speed [m/s]

# Returns
- `Diff_ra`: Residual in aerodynamic resistance [s/m]
"""
function solve_obhukov_length(
    LAN::FT, ra::FT, zom::FT, zoh::FT, z::FT, u::FT
) where {FT<:AbstractFloat}
    k = FT(0.4)  # von Karman constant

    # Calculate stability functions
    Fih_z, Fim_z = businger_stability_functions(z/LAN)
    _, Fim_zom = businger_stability_functions(zom/LAN)
    Fih_zoh, _ = businger_stability_functions(zoh/LAN)

    # Garrat 1992 - Aerodynamic resistance for heat flux [s/m]
    raMO = (1/(u * k^2)) * ((log(z/zom) - Fim_z + Fim_zom) * (log(z/zoh) - Fih_z + Fih_zoh))

    return raMO - ra
end
