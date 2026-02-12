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

    f(x, p) = solve_obhukov_length(x, ra, zom, zoh, z, u)
    rb = -eps(FT)
    frb = f(rb, nothing)
    lb = -10
    flb = f(lb, nothing)
    while sign(frb) == sign(flb)
        rb *= 10
        frb = f(rb, nothing)
        if abs(rb) > 1e6
            @debug "backcalculate_obhukov_length does not change sign in the range, setting LAN to -Inf."
            return FT(-Inf), solve_obhukov_length(FT(-Inf), ra, zom, zoh, z, u)
        end
    end

    prob = IntervalNonlinearProblem(
        (x, p) -> solve_obhukov_length(x, ra, zom, zoh, z, u), (lb, rb)
    )
    sol = solve(prob, Brent(); abstol=FT(1e-6), maxiters=400)

    if !successful_retcode(sol)
        throw(
            ErrorException(
                "Failed to converge when backcalculating Obhukov length. Solution return code: $(sol.retcode)",
            ),
        )
    end

    LAN = sol.u

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
function businger_stability_functions(y)
    if y > 0
        # Stable condition
        a = 1
        b = 0.667
        c = 5
        d = 0.35

        Fim = -(a * y + b * (y - c/d) * exp(-d * y) + b * c/d)
        Fih = -(((1 + 2 * a * y/3)^1.5) + b * (y - c/d) * exp(-d * y) + (b * c/d - 1))
    else
        # Unstable condition
        G = 16  # Dyer (1974)
        x = (1 - G * y)^(0.25)
        Fim = log((0.5 * (1 + x^2)) * ((0.5 * (1 + x))^2)) - 2 * atan(x) + π/2
        Fih = 2 * log(0.5 * (1 + x^2))
    end
    return Fih, Fim
end
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
    LAN, ra::FT, zom::FT, zoh::FT, z::FT, u::FT
) where {FT<:AbstractFloat}
    k = FT(0.4)  # von Karman constant
    # @info LAN

    # Calculate stability functions
    Fih_z, Fim_z = businger_stability_functions(z/LAN)
    _, Fim_zom = businger_stability_functions(zom/LAN)
    Fih_zoh, _ = businger_stability_functions(zoh/LAN)

    # Garrat 1992 - Aerodynamic resistance for heat flux [s/m]
    raMO = (1/(u * k^2)) * ((log(z/zom) - Fim_z + Fim_zom) * (log(z/zoh) - Fih_z + Fih_zoh))

    return raMO - ra
end
