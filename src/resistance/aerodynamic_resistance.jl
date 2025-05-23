"""
    aerodynamic_resistance(Ta, Ts, Pre, zatm, disp_h, zom, zoh, Ws, ea)

Calculate aerodynamic resistance from individual meteorological parameters.

# Arguments
- `Ta`: Air temperature [°C]
- `Ts`: Surface temperature [°C]
- `Pre`: Atmospheric pressure [hPa]
- `zatm`: Reference height [m]
- `disp_h`: Zero plane displacement height [m]
- `zom`: Roughness length for momentum [m]
- `zoh`: Roughness length for heat [m]
- `Ws`: Wind speed [m/s]
- `ea`: Vapor pressure [Pa]

# Returns
- `ra`: Aerodynamic resistance [s/m]
"""
function aerodynamic_resistance(
    Ta::FT, Ts::FT, Pre::FT, zatm::FT, disp_h::FT, zom::FT, zoh::FT, Ws::FT, ea::FT, es::FT
) where {FT<:AbstractFloat}
    if disp_h > zatm
        throw(ArgumentError("disp_h must be less zatm."))
    end

    # Convert pressure to Pa
    Pre *= FT(100)

    # Wind speed and measurement height
    u = Ws
    z = zatm - disp_h

    # Physical constants
    g = FT(9.81)                                    # Gravity acceleration [m/s²]
    P_Ref = FT(100000)                              # Reference pressure [Pa]
    Rd = FT(287.05)                                 # Dry air gas constant [J/kg/K]
    cp = 1005 + ((Ta + FT(23.15))^2) / FT(3364)     # Specific heat of air [J/kg/K]
    k = FT(0.4)                                     # Von Karman constant

    Oa = (Ta + FT(273.15)) * (Pre / P_Ref)^(-Rd / cp)  # Air potential temperature [K]
    Os = (Ts + FT(273.15)) * (Pre / P_Ref)^(-Rd / cp)  # Surface potential temperature [K]

    alpha = 1.0

    Ova = Oa * Pre ./ (Pre - FT(0.378) * ea)  # Air virtual potential temperature [K]
    Ovs = Os * Pre ./ (Pre - FT(0.378) * alpha * es)  # Surface virtual potential temperature [K]

    # Richardson number calculation
    f2 = ((1 - zom / z)^2) / (1 - zoh / z)  # Kot and Song (1998) adjustment
    Ri = -f2 * ((g * z * (Ovs - Ova)) / (FT(0.5) * (Ovs + Ova) * u^2))

    # Use approximate solution CH = CDN*f(Ri)
    # Following Mascart et al., 1995 and others
    mu = log(zom / zoh)
    ph = FT(0.5802) - FT(0.1571) * mu + FT(0.0327) * mu^2 - FT(0.0026) * mu^3
    Chs = FT(3.2165) + FT(4.3431) * mu + FT(0.5360) * mu^2 - FT(0.0781) * mu^3

    # Transport coefficients
    Cdn = (k^2) / (log(z / zom)^2)  # Neutral condition
    Ch = FT(15) * Chs * Cdn * (z / zoh)^ph * (log(z / zom) / log(z / zoh))

    # Stability correction
    if Ri <= zero(FT)
        # Unstable
        FH = (1 - FT(15) * Ri / (1 + Ch * sqrt(abs(Ri)))) * (log(z / zom) / log(z / zoh))
    else
        # Stable
        FH = (1 / (1 + FT(15) * Ri * sqrt(1 + FT(5) * Ri))) * (log(z / zom) / log(z / zoh))
    end

    # Final transport coefficient and resistance
    CH = FH * Cdn
    ra = 1 / (CH * u)

    # Handle free convection case (Beljaars 1994)
    if (u <= FT(0.05)) && (Ovs > Ova)
        Cs = FT(0.15)
        ni = FT(1.51e-5)  # Air kinematic viscosity [m²/s]
        Pr = FT(0.71)     # Prandtl number
        CHu =
            Cs *
            (g * ni / (FT(0.5) * (Ovs + Ova) * Pr^2))^(FT(1) / 3) *
            (Ovs - Ova)^(FT(1) / 3)
        ra = 1 / CHu
    end

    return ra
end
