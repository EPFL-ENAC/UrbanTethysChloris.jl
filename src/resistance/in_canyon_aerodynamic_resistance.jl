"""
    in_canyon_aerodynamic_resistance(
        uatm::FT,
        Zatm::FT,
        Ta::FT,
        Ts::FT,
        hcan_max::FT,
        hcan::FT,
        dcan::FT,
        zomcan::FT,
        Zref_und::FT,
        zom_und::FT,
        Zp1::FT,
        Zp2::FT,
        Zp3::FT,
        Pre::FT,
        ea::FT,
        RoughnessParameter::Symbol
    ) where {FT<:AbstractFloat}

Calculate aerodynamic resistance within urban canyon.

# Arguments
- `uatm`: Wind speed at atmospheric reference height [m/s]
- `Zatm`: Atmospheric reference height [m]
- `Ta`: Air temperature [°C]
- `Ts`: Surface temperature [°C]
- `hcan_max`: Maximum canyon height [m]
- `hcan`: Urban canyon height [m]
- `dcan`: Urban canyon displacement height [m]
- `zomcan`: Urban roughness length [m]
- `Zref_und`: Reference height within canyon [m]
- `zom_und`: Roughness length of ground surface [m]
- `Zp1`: First height within canyon [m]
- `Zp2`: Second height within canyon [m]
- `Zp3`: Third height within canyon [m]
- `Pre`: Atmospheric pressure [Pa]
- `ea`: Vapor pressure [Pa]
- `RoughnessParameter`: Roughness parameterization type

# Returns
- `rap_can`: Urban undercanopy resistance from ground to displacement height [s/m]
- `rap_Zp1`: Undercanopy resistance from ground to Zp1 [s/m]
- `rap_Zp1_In`: Inverse resistance from Zp1 to displacement height [s/m]
- `rap_Zp2`: Undercanopy resistance from ground to Zp2 [s/m]
- `rap_Zp2_In`: Inverse resistance from Zp2 to displacement height [s/m]
- `rap_Zp3`: Undercanopy resistance from ground to Zp3 [s/m]
- `rap_Zp3_In`: Inverse resistance from Zp3 to displacement height [s/m]
- `u_Hcan`: Wind speed at canyon height [m/s]
- `u_Zp1`: Wind speed at Zp1 [m/s]
- `u_Zp2`: Wind speed at Zp2 [m/s]
- `u_Zp3`: Wind speed at Zp3 [m/s]
- `uref_und`: Wind speed at reference height [m/s]
- `alpha`: Canyon attenuation coefficient [-]
"""
function in_canyon_aerodynamic_resistance(
    uatm::FT,
    Zatm::FT,
    Ta::FT,
    Ts::FT,
    hcan_max::FT,
    hcan::FT,
    dcan::FT,
    zomcan::FT,
    Zref_und::FT,
    zom_und::FT,
    Zp1::FT,
    Zp2::FT,
    Zp3::FT,
    Pre::FT,
    ea::FT,
    RoughnessParameter::Symbol,
) where {FT<:AbstractFloat}
    # Adjust reference heights if needed
    if Zref_und > (dcan + zomcan)
        @info "The urban undercanopy reference height is too big or the canyon too shallow: Zref_und>(dcan+zomcan)"
        Zref_und = dcan + zomcan - FT(0.01)
    end
    if zom_und > Zref_und
        @info "The ground surface roughness is bigger than the urban undercanopy reference height: zom_und>Zref_und"
        zom_und = Zref_und - FT(0.01)
    end
    if RoughnessParameter == :MacD
        hcan_max = hcan
    end

    # Adjust calculation heights
    Zp1 = clamp(Zp1, Zref_und, dcan + zomcan)
    Zp2 = clamp(Zp2, Zref_und, dcan + zomcan)
    Zp3 = clamp(Zp3, Zref_und, dcan + zomcan)

    # Constants
    g = FT(9.81)    # Gravity [m/s²]
    k = FT(0.41)    # von Karman constant

    # Wind speeds calculations
    us_atm = k * uatm / log((Zatm - dcan) / zomcan)  # Friction velocity [m/s]
    u_Hcan = (us_atm / k) * log((hcan_max - dcan) / zomcan)  # Wind speed at canyon top [m/s]
    alpha = log(uatm / u_Hcan) / (Zatm / hcan_max - 1)  # Attenuation coefficient
    uref_und = u_Hcan * exp(-alpha * (1 - Zref_und / hcan_max))  # Reference wind speed
    Kh_can = k^2 * uatm * (hcan_max - dcan) / log((Zatm - dcan) / zomcan)  # Eddy diffusion

    # Canyon resistance calculation
    rap_can =
        hcan_max * exp(alpha) / (Kh_can * alpha) *
        (exp(-alpha * (Zref_und / hcan_max)) - exp(-alpha * ((dcan + zomcan) / hcan_max))) +
        1 / (k^2 * uref_und) * log(Zref_und / zom_und)^2

    # Stability correction
    P_Ref = FT(100000)  # Reference pressure [Pa]
    Rd = FT(287.05)     # Gas constant [J/kg/K]
    cp = FT(1005) + ((Ta + FT(23.15))^2) / FT(3364)  # Specific heat [J/kg/K]

    # Potential temperatures
    Oa = (Ta + FT(273.15)) * (Pre/P_Ref)^(-Rd/cp)  # Air [K]
    Os = (Ts + FT(273.15)) * (Pre/P_Ref)^(-Rd/cp)  # Surface [K]
    Ova = Oa
    Ovs = Os

    # Richardson number calculation
    Ri2 = (g * (Ova - Ovs) * Zref_und) / (uref_und^2 * (FT(0.5) * (Ova + Ovs)))
    Ri2 = min(Ri2, FT(0.16))  # Max stability

    # Apply stability correction
    rap_can *= if Ri2 < 0
        1 / (1 - 5 * Ri2)^(FT(0.75))  # Unstable
    else
        1 / (1 - 5 * Ri2)^2  # Stable
    end

    # Calculate for each height level
    function process_level(Zp)
        u_Zp = u_Hcan * exp(-alpha * (1 - Zp/hcan_max))
        rap_Zp =
            hcan_max * exp(alpha) / (Kh_can * alpha) *
            (exp(-alpha * (Zref_und/hcan_max)) - exp(-alpha * (Zp/hcan_max))) +
            1 / (k^2 * uref_und) * log(Zref_und/zom_und)^2

        rap_Zp *= if Ri2 < 0
            1 / (1 - 5 * Ri2)^(FT(0.75))
        else
            1 / (1 - 5 * Ri2)^2
        end

        rap_Zp_In = max(rap_can - rap_Zp, FT(0.1))

        return u_Zp, rap_Zp, rap_Zp_In
    end

    # Process each level
    u_Zp1, rap_Zp1, rap_Zp1_In = process_level(Zp1)
    u_Zp2, rap_Zp2, rap_Zp2_In = process_level(Zp2)
    u_Zp3, rap_Zp3, rap_Zp3_In = process_level(Zp3)

    return rap_can,
    rap_Zp1,
    rap_Zp1_In,
    rap_Zp2,
    rap_Zp2_In,
    rap_Zp3,
    rap_Zp3_In,
    u_Hcan,
    u_Zp1,
    u_Zp2,
    u_Zp3,
    uref_und,
    alpha
end
