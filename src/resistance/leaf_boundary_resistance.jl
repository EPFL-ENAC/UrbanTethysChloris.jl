"""
    leaf_boundary_resistance(u_hc, Ts, Ta, d_leaf, alpha)
    leaf_boundary_resistance(Ws, Ts, Ta, hc, d_leaf, zatm, disp_h, zom)

Calculate leaf boundary layer resistance following Vesala (1998) and Choudhury & Monteith (1988).

# Arguments
- `u_hc`: Wind speed at canopy top [m/s]
- `Ts`: Surface temperature [°C]
- `Ta`: Air temperature [°C]
- `d_leaf`: Leaf dimension [cm]
- `alpha`: Attenuation coefficient [-]
- `Ws`: Wind speed at height [m/s]
- `hc`: Canopy height [m]
- `zatm`: Atmospheric height [m]
- `disp_h`: Displacement height [m]
- `zom`: Roughness length [m]

# Returns
- `rb`: Leaf boundary layer resistance [s/m] for one-sided unit leaf area
"""

function leaf_boundary_resistance(u_hc, Ts, Ta, d_leaf, alpha)
    # Convert leaf dimension to meters
    d_leaf = d_leaf / 100  # [m]

    # Constants
    a = 0.005    # [m/s^0.5] --  Chodhury and Monteith 1988
    Dh = 1.9e-5  # [m2/s]

    # Calculate forced convection component
    gb = (2 * a / alpha) * sqrt(u_hc / d_leaf) * (1 - exp(-alpha / 2))  # [m/s]

    # Calculate free convection component
    Gr = 1.6e8 * (Ts - Ta) * (d_leaf^3) * (Ts > Ta)  # [-]
    gb_free = 0.5 * Dh * (Gr^0.25) / d_leaf      # [m/s]

    # Combine and convert to resistance
    gb = gb + gb_free
    rb = 1 / gb  # Leaf Boundary Layer Resistance [s/m] one-sided for unit leaf

    return rb
end

function leaf_boundary_resistance(Ws, Ts, Ta, hc, d_leaf, LAI, zatm, disp_h, zom)
    u = Ws
    k = 0.4
    d = disp_h
    z = zatm

    # Hypothesis Logarithmic distribution of wind speed
    us = k * u / log((z - d) / zom)
    u_hc = (us / k) * log((hc - d) / zom)

    # Attenuation coefficient
    alpha = log(u / u_hc) / (z / hc - 1)
    alpha *= 0.5 * LAI / 2

    return leaf_boundary_resistance(u_hc, Ts, Ta, d_leaf, alpha)
end
