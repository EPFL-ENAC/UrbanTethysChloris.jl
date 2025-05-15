"""
    soil_parameters(Psan::FT, Pcla::FT, Porg::FT) where {FT<:AbstractFloat}

Calculate soil hydraulic and thermal parameters based on soil texture composition.

# Arguments
- `Psan::AbstractFloat`: Sand fraction in soil [0-1]
- `Pcla::AbstractFloat`: Clay fraction in soil [0-1]
- `Porg::AbstractFloat`: Organic matter fraction in soil [0-1]

# Returns
A tuple containing:
- `Osat`: Saturated soil water content [-]
- `L`: Slope of logarithmic tension-moisture curve [-]
- `Pe`: Air entry tension [kPa]
- `Ks`: Saturated hydraulic conductivity [mm/h]
- `O33`: Soil water content at 33 kPa [-]
- `rsd`: Soil bulk density [kg/m³]
- `lan_dry`: Thermal conductivity of dry soil [W/m K]
- `lan_s`: Thermal conductivity of solid soil components [W/m K]
- `cv_s`: Volumetric heat capacity of solid soil components [J/m³ K]
- `K_usle`: USLE K factor for soil erodibility [kg*h/J*mm]
"""
function soil_parameters(Psan::FT, Pcla::FT, Porg::FT) where {FT<:AbstractFloat}
    # Check if the fractions are within the valid range
    if max(Psan, Pcla, Porg) > 1 || min(Psan, Pcla, Porg) < 0
        throw(ArgumentError("Fractions must be between 0 and 1"))
    end

    # Initialize variables
    Rw = 0.0    # Weight fraction of gravel [g Gravel /g bulk soil]
    DF = 1.0    # density factor

    Psil = 1.0 - Psan - Pcla - Porg
    if Psil < 0
        throw(ArgumentError("Sum of fractions must be ≤ 1"))
    end

    Porg = Porg * 100

    # Parameterization [Saxton and Rawls 2006]
    O1500t =
        -0.024 * Psan + 0.487 * Pcla + 0.006 * Porg + 0.005 * (Psan * Porg) -
        0.013 * (Pcla * Porg) +
        0.068 * (Psan * Pcla) +
        0.031
    O1500 = O1500t + 0.14 * O1500t - 0.02    # 1500 kPa Moisture

    O33t =
        -0.251 * Psan + 0.195 * Pcla + 0.011 * Porg + 0.006 * (Psan * Porg) -
        0.027 * (Pcla * Porg) +
        0.452 * (Psan * Pcla) +
        0.299
    O33 = O33t + (1.283 * (O33t)^2 - 0.374 * (O33t) - 0.015)    # 33 kPa Moisture

    Os_33t =
        0.278 * Psan + 0.034 * Pcla + 0.022 * Porg - 0.018 * (Psan * Porg) -
        0.027 * (Pcla * Porg) - 0.584 * (Psan * Pcla) + 0.078
    Os_33 = Os_33t + (0.636 * Os_33t - 0.107)    # SAT-33 kPa Moisture

    # Parameters
    B = (log(1500) - log(33)) / (log(O33) - log(O1500))    # Coefficient of moisture tension
    L = 1 / B    # Slope of logaritimc tension-moisture curve

    Osat = O33 + Os_33 - 0.097 * Psan + 0.043    # Saturation moisture 0 kPa
    rsd = (1 - Osat) * 2650    # normal density dry soil  [kg/m^3]

    # Density effects
    rsd_df = rsd * DF    # adj density density dry soil [kg/m^3]
    Osat_df = 1 - (rsd_df / 2650)    # Saturation moisture 0 kPa
    O33_df = O33 - 0.2 * (Osat - Osat_df)    # 33 kPa Moisture
    Os_33_df = Osat_df - O33_df    # SAT-33 kPa Moisture

    Osat = Osat_df
    O33 = O33_df
    Os_33 = Os_33_df
    rsd = rsd_df

    Ks = 1930 * (Osat - O33)^(3 - L)    # saturation conductivty [mm/h]
    Pet =
        -21.67 * Psan - 27.93 * Pcla - 81.97 * Os_33 +
        71.12 * (Psan * Os_33) +
        8.29 * (Pcla * Os_33) +
        14.05 * (Psan * Pcla) +
        27.16
    Pe = Pet + (0.02 * Pet^2 - 0.113 * Pet - 0.70)    # Tension at air entry [kPa]
    Pe = max.(Pe, 0.5)

    # Gravel Effects
    alpha = rsd / 2650    # matric soil density / gravel density
    Rv = (alpha * Rw) / (1 - Rw * (1 - alpha))    # Volume fraction of gravel [g/cm^3]
    rhoB = (rsd / 1000) * (1 - Rv) + (Rv * 2.65)    # Bulk soil density [g/cm^3]
    Osatb = Osat * (1 - Rv)
    O33b = O33 * (1 - Rv)
    Kb = Ks * (1 - Rw) / (1 - Rw * (1 - 3 * alpha / 2))    # Saturated conductivity bulk soil [mm/h]

    rsd = rhoB * 1000    # [kg/m^3]
    Osat = Osatb    # [-]
    O33 = O33b    # [-]
    Ks = Kb    # [mm/h]

    # Check consistency
    if !(Osat ≥ O33 ≥ O1500)
        error("Soil parameters are inconsistent")
    end

    # Thermal characteristics
    Osat_min = 0.489 - 0.1267 * Psan
    fom = (Osat - Osat_min) / (0.9 - Osat_min)
    fom = clamp(fom, 0, 1)

    lan_s_om = 0.25    # Thermal conductivity Organic matter [W/m K]
    lan_dry_om = 0.05    # Thermal conductivity Organic matter dry [W/m K]
    cv_s_om = 2.5e6    # Volumetric heat capacity Organic matter [J/m^3 K]

    rsd_s = 2700 * (1 - Osat)
    lan_dry_min = (0.135 * rsd_s + 64.7) / (2700 - 0.947 * rsd_s)
    lan_s_min = (8.8 * Psan + 2.92 * Pcla) / (Psan + Pcla)

    lan_dry = (1 - fom) * lan_dry_min + fom * lan_dry_om
    lan_s = (1 - fom) * lan_s_min + fom * lan_s_om
    cv_s_min = 1e6 * (2.128 * Psan + 2.385 * Pcla) / (Psan + Pcla)
    cv_s = (1 - fom) * cv_s_min + fom * cv_s_om

    # K USLE Parameter [Williams 1995]
    Porg_c = Porg / 1.72
    fsand = 0.2 + 0.3 * exp(-25.6 * Psan * (1 - Psil))
    fcli = (Psil / (Pcla + Psil))^0.3
    forg = 1 - (0.25 * Porg_c) / (Porg_c + exp(3.72 - 2.95 * Porg_c))
    fhisand = 1 - (0.7 * (1 - Psan)) / ((1 - Psan) + exp(-5.51 - 22.9 * (1 - Psan)))
    K_usle = fsand * fcli * forg * fhisand / 1000    # [kg*h/J*mm] erosivity factor

    return Osat, L, Pe, Ks, O33, rsd, lan_dry, lan_s, cv_s, K_usle
end
