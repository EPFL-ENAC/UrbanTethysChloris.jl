"""
    soil_parameters(
        Psan::FT,
        Pcla::FT,
        Porg::FT
    ) where {FT<:AbstractFloat}

Calculate soil hydraulic and thermal parameters based on soil texture composition.

# Arguments
- `Psan`: Sand fraction in soil [0-1]
- `Pcla`: Clay fraction in soil [0-1]
- `Porg`: Organic matter fraction in soil [0-1]

# Returns
- `Osat::FT`: Saturated soil water content [-]
- `L::FT`: Slope of logarithmic tension-moisture curve [-]
- `Pe::FT`: Air entry tension [kPa]
- `Ks::FT`: Saturated hydraulic conductivity [mm/h]
- `O33::FT`: Soil water content at 33 kPa [-]
- `rsd::FT`: Soil bulk density [kg/m³]
- `lan_dry::FT`: Thermal conductivity of dry soil [W/m K]
- `lan_s::FT`: Thermal conductivity of solid soil components [W/m K]
- `cv_s::FT`: Volumetric heat capacity of solid soil components [J/m³ K]
- `K_usle::FT`: USLE K factor for soil erodibility [kg*h/J*mm]

"""
function soil_parameters(Psan::FT, Pcla::FT, Porg::FT) where {FT<:AbstractFloat}
    # Check if the fractions are within the valid range
    if max(Psan, Pcla, Porg) > 1||min(Psan, Pcla, Porg) < 0
        throw(ArgumentError("Fractions must be between 0 and 1"))
    end

    # Initialize variables
    Rw = 0    # Weight fraction of gravel [g Gravel /g bulk soil]
    DF = FT(1)   # density factor

    Psil = 1 - Psan - Pcla - Porg
    if Psil < 0
        throw(ArgumentError("Sum of fractions must be ≤ 1"))
    end

    Porg = Porg * 100

    # Parameterization [Saxton and Rawls 2006]
    O1500t =
        FT(-0.024) * Psan +
        FT(0.487) * Pcla +
        FT(0.006) * Porg +
        FT(0.005) * (Psan * Porg) - FT(0.013) * (Pcla * Porg) +
        FT(0.068) * (Psan * Pcla) +
        FT(0.031)
    O1500 = O1500t + FT(0.14) * O1500t - FT(0.02)    # 1500 kPa Moisture
    O33t =
        FT(-0.251) * Psan +
        FT(0.195) * Pcla +
        FT(0.011) * Porg +
        FT(0.006) * (Psan * Porg) - FT(0.027) * (Pcla * Porg) +
        FT(0.452) * (Psan * Pcla) +
        FT(0.299)
    O33 = O33t + (FT(1.283) * (O33t)^2 - FT(0.374) * (O33t) - FT(0.015))    # 33 kPa Moisture

    Os_33t =
        FT(0.278) * Psan + FT(0.034) * Pcla + FT(0.022) * Porg - FT(0.018) * (Psan * Porg) -
        FT(0.027) * (Pcla * Porg) - FT(0.584) * (Psan * Pcla) + FT(0.078)
    Os_33 = Os_33t + (FT(0.636) * Os_33t - FT(0.107))    # SAT-33 kPa Moisture

    # Parameters
    B = (log(FT(1500)) - log(FT(33))) / (log(O33) - log(O1500))    # Coefficient of moisture tension
    L = 1 / B    # Slope of logaritimc tension-moisture curve

    Osat = O33 + Os_33 - FT(0.097) * Psan + FT(0.043)    # Saturation moisture 0 kPa
    rsd = (1 - Osat) * 2650    # normal density dry soil  [kg/m^3]

    # Density effects
    rsd_df = rsd * DF    # adj density density dry soil [kg/m^3]
    Osat_df = 1 - (rsd_df / 2650)    # Saturation moisture 0 kPa
    O33_df = O33 - FT(0.2) * (Osat - Osat_df)    # 33 kPa Moisture
    Os_33_df = Osat_df - O33_df    # SAT-33 kPa Moisture

    Osat = Osat_df
    O33 = O33_df
    Os_33 = Os_33_df
    rsd = rsd_df

    Ks = 1930 * (Osat - O33)^(3 - L)    # saturation conductivty [mm/h]
    Pet =
        FT(-21.67) * Psan - FT(27.93) * Pcla - FT(81.97) * Os_33 +
        FT(71.12) * (Psan * Os_33) +
        FT(8.29) * (Pcla * Os_33) +
        FT(14.05) * (Psan * Pcla) +
        FT(27.16)
    Pe = Pet + (FT(0.02) * Pet^2 - FT(0.113) * Pet - FT(0.70))    # Tension at air entry [kPa]
    Pe = max.(Pe, FT(0.5))

    # Gravel Effects
    alpha = rsd / 2650   # matric soil density / gravel density
    Rv = (alpha * Rw) / (1 - Rw * (1 - alpha))    # Volume fraction of gravel [g/cm^3]
    rhoB = (rsd / 1000) * (1 - Rv) + (Rv * FT(2.65))    # Bulk soil density [g/cm^3]
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

    # Thermal characteristics of soil [Farouki 1981] [Williams 1963] [Oleson et al. 2004]
    rsd_s = 2700 * (1 - Osat)
    lan_dry = (FT(0.135) * rsd_s + FT(64.7)) / (2700 - FT(0.947) * rsd_s)
    lan_s = (FT(8.8) * Psan + FT(2.92) * Pcla) / (Psan + Pcla)
    cv_s = FT(1e6) * (FT(2.128) * Psan + FT(2.385) * Pcla) / (Psan + Pcla)

    # K USLE Parameter [Williams 1995]
    Porg_c = Porg / FT(1.72)
    fsand = FT(0.2) + FT(0.3) * exp(FT(-25.6) * Psan * (1 - Psil))
    fcli = (Psil / (Pcla + Psil))^FT(0.3)
    forg = 1 - (FT(0.25) * Porg_c) / (Porg_c + exp(FT(3.72) - FT(2.95) * Porg_c))
    fhisand =
        1 - (FT(0.7) * (1 - Psan)) / ((1 - Psan) + exp(FT(-5.51) - FT(22.9) * (1 - Psan)))
    K_usle = fsand * fcli * forg * fhisand / 1000    # [kg*h/J*mm] erosivity factor

    return Osat, L, Pe, Ks, O33, rsd, lan_dry, lan_s, cv_s, K_usle
end
