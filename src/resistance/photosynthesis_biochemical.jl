
"""
    photosynthesis_biochemical(Cc, IPAR, Csl, ra, rb, Ts, Pre, Ds, Psi_L, Psi_sto_50,
                             Psi_sto_00, CT, Vmax, DSE, Ha, FI, Oa, Do, a1, go, gmes, rjv)

Calculate photosynthesis biochemical processes at leaf level.

# Arguments
- `Cc`: CO₂ concentration in mesophyll cells [μmolCO₂/molAIR]
- `IPAR`: Intercepted photosynthetically active radiation [W/m²]
- `Csl`: Leaf surface CO₂ concentration [μmolCO₂/molAIR]
- `ra`: Aerodynamic resistance [s/m]
- `rb`: Boundary layer resistance [s/m]
- `Ts`: Leaf temperature [°C]
- `Pre`: Atmospheric pressure [hPa]
- `Ds`: Vapor pressure deficit
- `Psi_L`: Leaf water potential
- `Psi_sto_50`: Water potential at 50% stomatal closure
- `Psi_sto_00`: Water potential at complete stomatal closure
- `CT`: Photosynthesis type (3=C3, 4=C4)
- `Vmax`: Maximum Rubisco capacity at 25°C [μmolCO₂/s/m²]
- `DSE`: Entropy factor [kJ/mol/K]
- `Ha`: Activation energy [kJ/mol]
- `FI`: Intrinsic quantum efficiency [μmolCO₂/μmolPhotons]
- `Oa`: Intercellular O₂ partial pressure [μmolO₂/mol]
- `Do`: Empirical coefficient for vapor pressure
- `a1`: Empirical coefficient for stomatal conductance
- `go`: Minimum stomatal conductance [molCO₂/s/m²]
- `gmes`: Mesophyll conductance [molCO₂/s/m²]
- `rjv`: Ratio of Jmax to Vmax

# Returns
- `CcF`: Final CO₂ concentration in mesophyll [μmolCO₂/molAIR]
- `An`: Net assimilation rate [μmolCO₂/s/m²]
- `rs`: Stomatal resistance [s/m]
- `Rdark`: Dark respiration [μmolCO₂/s/m²]
- `F755nm`: Chlorophyll fluorescence [W/m²/sr/μm]
- `GAM`: CO₂ compensation point [Pa]
- `gsCO2`: Stomatal conductance to CO₂ [μmolCO₂/s/m²]
"""
function photosynthesis_biochemical(
    Cc::FT,
    IPAR::FT,
    Csl::FT,
    ra::FT,
    rb::FT,
    Ts::FT,
    Pre::FT,
    Ds::FT,
    Psi_L::FT,
    Psi_sto_50::FT,
    Psi_sto_00::FT,
    CT::Int,
    Vmax::FT,
    DSE::FT,
    Ha::FT,
    FI::FT,
    Oa::FT,
    Do::FT,
    a1::FT,
    go::FT,
    gmes::FT,
    rjv::FT,
) where {FT<:AbstractFloat}
    Ta = Ts
    Pre0 = FT(101325)  # Reference pressure [Pa]
    Tf = FT(273.15)    # Freezing point [K]

    # Unit conversions
    Pre *= 100  # [Pa]
    IPAR *= FT(4.57)  # [μmolPhotons/s/m²]
    ra *= (FT(0.0224) * (Ta + Tf) * Pre0 / (Tf * Pre)) * FT(1e-6)  # [m²s/μmolH₂O]
    rb *= (FT(0.0224) * (Ta + Tf) * Pre0 / (Tf * Pre)) * FT(1e-6)  # [m²s/μmolH₂O]
    Cc *= FT(1e-6) * Pre   # [Pa]
    Oa *= FT(1e-6) * Pre   # [Pa]
    Csl *= FT(1e-6) * Pre  # [Pa]

    # Convert conductances
    rmes = FT(1) / (FT(1e6) * gmes)  # [s m²/μmolCO₂]
    go *= FT(1e6)  # [μmolCO₂/s/m²]

    # Temperature dependencies
    Ts_k = Ts + Tf  # [K]
    Tref = FT(25) + Tf  # [K]
    R = FT(0.008314)  # Gas constant [kJ/K/mol]

    # Kattge and Knorr 2007 temperature dependencies
    Hd = FT(200)  # Deactivation energy [kJ/mol]
    kT =
        exp(Ha * (Ts_k - Tref) / (Tref * R * Ts_k)) *
        (1 + exp((Tref * DSE - Hd) / (Tref * R))) /
        (1 + exp((Ts_k * DSE - Hd) / (Ts_k * R)))
    Vm = Vmax * kT

    # Maximum electron transport rate
    Jmax = Vmax * rjv
    Ha_J = FT(50)
    DS_J = FT(0.646)
    kT_J =
        exp(Ha_J * (Ts_k - Tref) / (Tref * R * Ts_k)) *
        (1 + exp((Tref * DS_J - Hd) / (Tref * R))) /
        (1 + exp((Ts_k * DS_J - Hd) / (Ts_k * R)))
    Jm = Jmax * kT_J

    # Triose phosphate utilization
    Ha_T = FT(53.1)
    DS_T = FT(0.490)
    Hd_T = FT(150.65)
    TPU25 = FT(0.1182) * Vmax
    kT_T =
        exp(Ha_T * (Ts_k - Tref) / (Tref * R * Ts_k)) *
        (1 + exp((Tref * DS_T - Hd_T) / (Tref * R))) /
        (1 + exp((Ts_k * DS_T - Hd_T) / (Ts_k * R)))
    TPU = TPU25 * kT_T

    # C4 specific temperature dependencies
    if CT == 4
        s1 = FT(0.3)
        s3 = FT(0.2)
        Tup = FT(40)
        Tlow = FT(15)
        f1T = 1 / (1 + exp(s1 * (Ts - Tup)))
        f2T = 1 / (1 + exp(s3 * (Tlow - Ts)))
        fT = FT(2.0)^(FT(0.1) * (Ts - 25))
        Vm = Vmax * fT * f1T * f2T
        ke25 = FT(20000) * Vmax
        ke = ke25 * fT
    end

    # CO₂ compensation point (Bonan et al 2011)
    Ha_G = FT(37.83)
    kT_G = exp(Ha_G * (Ts_k - Tref) / (Tref * R * Ts_k))
    GAM25 = FT(42.75) * FT(1e-6) * Pre
    GAM = GAM25 * kT_G

    # Michaelis-Menten constants for C3
    if CT == 3
        Ha_Kc = FT(79.43)
        Kc25 = FT(404.9) * FT(1e-6) * Pre
        kT_Kc = exp(Ha_Kc * (Ts_k - Tref) / (Tref * R * Ts_k))
        Kc = Kc25 * kT_Kc

        Ha_Ko = FT(36.38)
        Ko25 = FT(278.4) * FT(1e-3) * Pre
        kT_Ko = exp(Ha_Ko * (Ts_k - Tref) / (Tref * R * Ts_k))
        Ko = Ko25 * kT_Ko
    end

    # Dark respiration
    if CT == 3
        # Bonan et al 2011
        Ha_Rd = FT(46.39)
        DS_Rd = FT(0.490)
        Hd_Rd = FT(150.65)
        Rdark25 = FT(0.015) * Vmax
        # kT_Rd = exp(Ha_Rd*(Ts_k-Tref)/(Tref*R*Ts_k))
        kT_Rd =
            exp(Ha_Rd * (Ts_k - Tref) / (Tref * R * Ts_k)) *
            (1 + exp((Tref * DS_Rd - Hd_Rd) / (Tref * R))) /
            (1 + exp((Ts_k * DS_Rd - Hd_Rd) / (Ts_k * R)))
        Rdark = Rdark25 * kT_Rd

    elseif CT == 4
        # Collatz 1991 - 1992 - Sellers 1996 Bonan et al 2011
        fT = FT(2.0)^(FT(0.1) * (Ts - 25))
        fT3 = 1 / (1 + exp(1.3 * (Ts - 55)))
        Rdark25 = FT(0.025) * Vmax
        Rdark = Rdark25 * fT * fT3
    end

    # Photosynthesis factors
    Q = FI * IPAR  # Light absorbed by PSII in CO₂ units
    d1 = FT(0.7)   # d1 = 0.95 Leuning 1995; d1 = 0.7 Bonan et al., 2011
    d2 = -(Q + Jm / 4)
    d3 = Q * Jm / 4
    # Note: original MATLAB uses min(roots([d1,d2,d3]))
    J = min(
        (-d2 + sqrt(d2^2 - 4 * d1 * d3)) / (2 * d1),
        (-d2 - sqrt(d2^2 - 4 * d1 * d3)) / (2 * d1),
    )

    # Calculate assimilation rates
    if CT == 3
        JC = Vm * (Cc - GAM) / (Cc + Kc * (1 + Oa / Ko))
        JL = J * (Cc - GAM) / (4 * (Cc + FT(2) * GAM))
        JE = TPU * FT(3)
    elseif CT == 4
        JC = Vm
        JL = J / FT(4)
        JE = Tke * Cc / Pre
    end

    # First Polynomium
    b1 = CT == 3 ? FT(0.98) : FT(0.80)
    b2 = -(JC + JL)
    b3 = JC * JL
    JP = min(
        (-b2 + sqrt(b2^2 - 4 * b1 * b3)) / (2 * b1),
        (-b2 - sqrt(b2^2 - 4 * b1 * b3)) / (2 * b1),
    )

    # Second Polynomium
    c1 = FT(0.95)
    c2 = -(JP + JE)
    c3 = JP * JE
    A = min(
        (-c2 + sqrt(c2^2 - 4 * c1 * c3)) / (2 * c1),
        (-c2 - sqrt(c2^2 - 4 * c1 * c3)) / (2 * c1),
    )

    # Water stress function
    Rgsws = FT(0.02)
    p2 = log((1 - Rgsws) / Rgsws) / (Psi_sto_00 - Psi_sto_50) # [1/MPa]
    q2 = -p2 * Psi_sto_50  # [-]
    Rgsw = 1 / (1 + exp(p2 * Psi_L + q2))  # [fraction]
    fO = (1 - Rgsw)
    fO = clamp(fO, zero(FT), one(FT))

    # Solar-induced chlorophyll fluorescence (SIF)
    # Calculate electron transport from CO2 exchange
    Jfe = if CT == 3
        A * (Cc + 2 * GAM) / (Cc - GAM)
    else  # CT == 4
        A  # [μmolCO₂/s/m²]
    end

    # Fluorescence calculation
    fiP0 = FI * 4  # [μmol Electrons/μmolPhotons]
    fiP = fiP0 * Jfe / Q  # [0.4 max - stress decrease]
    dls = 1 - fiP / fiP0  # degree of light saturation

    kf = FT(0.05)
    kd = max(FT(0.03) * Ts + FT(0.0773), FT(0.087))
    kn = (FT(6.2473) * dls - FT(0.5944)) * dls

    fiF = kf / (kf + kd + kn) * (1 - fiP)  # [μmol Electrons/μmolPhotons]
    SIF = IPAR * fiF  # [μmol electrons/s m²]

    # Scale fluorescence
    k = FT(0.0375) * Vmax + FT(8.25)  # [μmol m⁻² s⁻¹ / W m⁻² sr⁻¹ μm⁻¹]
    F755nm = SIF / k  # [W m⁻² sr⁻¹ μm⁻¹]

    # Apply stress to assimilation
    A *= fO # Gross Assimilation Rate [umolCO2/ s m^2 ]

    # Net assimilation and conductance
    An = A - Rdark  # Net Assimilation Rate [μmolCO₂/s/m²]

    # Ball-Woodrow-Berry model with Dewar (2002) correction
    gsCO2 = max(go, go + a1 * An * Pre / ((Cc - GAM) * (1 + Ds / Do)))  # [μmolCO₂/s/m²]
    rsCO2 = 1 / gsCO2  # [s m²/μmolCO₂]

    # Calculate CO₂ concentration
    CcF = Csl - An * Pre * (rsCO2 + rmes + FT(1.37) * rb + ra)  # [Pa]
    CcF = max(CcF, zero(FT))

    # Calculate stomatal resistance
    rsH2O = (rsCO2 / FT(1.64)) * FT(1e6)  # [s/m²/molH₂O]
    An = (Csl - CcF) / (Pre * (rsCO2 + rmes + FT(1.37) * rb + ra))  # [μmolCO₂/s/m²]

    CcF /= (Pre * FT(1e-6))  # Convert back to [μmolCO₂/molAIR]
    rs = rsH2O * (Tf * Pre) / (FT(0.0224) * (Ts + Tf) * Pre0)  # [s/m]

    return CcF, An, rs, Rdark, F755nm, GAM, gsCO2
end
