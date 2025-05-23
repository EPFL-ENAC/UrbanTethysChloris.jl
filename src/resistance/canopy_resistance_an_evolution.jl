"""
    canopy_resistance_an_evolution(
        PAR_sun::FT,
        PAR_shd::FT,
        LAI::FT,
        Kopt::FT,
        Knit::FT,
        Fsun::FT,
        Fshd::FT,
        Citm1_sun::FT,
        Citm1_shd::FT,
        Ca::FT,
        ra::FT,
        rb::FT,
        Ts::FT,
        Ta::FT,
        Pre::FT,
        Ds::FT,
        Psi_L::FT,
        Psi_sto_50::FT,
        Psi_sto_99::FT,
        CT::Int,
        Vmax::FT,
        DS::FT,
        Ha::FT,
        FI::FT,
        Oa::FT,
        Do::FT,
        a1::FT,
        go::FT,
        e_rel::FT,
        e_relN::FT,
        gmes::FT,
        rjv::FT,
        mSl::FT,
        Sl::FT,
        atol::FT
    ) where {FT<:AbstractFloat}

Calculate canopy resistance and assimilation.

# Arguments
- `PAR_sun`: Photosynthetically active radiation for sunlit fraction [W/m²]
- `PAR_shd`: Photosynthetically active radiation for shaded fraction [W/m²]
- `LAI`: Leaf area index [-]
- `Kopt`: Optical extinction coefficient [-]
- `Knit`: Nitrogen extinction coefficient [-]
- `Fsun`: Sunlit fraction [-]
- `Fshd`: Shaded fraction [-]
- `Citm1_sun`: Previous CO₂ concentration for sunlit fraction [μmolCO₂/mol]
- `Citm1_shd`: Previous CO₂ concentration for shaded fraction [μmolCO₂/mol]
- `Ca`: Atmospheric CO₂ concentration [μmolCO₂/mol]
- `ra`: Aerodynamic resistance [s/m]
- `rb`: Boundary layer resistance [s/m]
- `Ts`: Surface temperature [°C]
- `Ta`: Air temperature [°C]
- `Pre`: Atmospheric pressure [hPa]
- `Ds`: Vapor pressure deficit [Pa]
- `Psi_L`: Leaf water potential [MPa]
- `Psi_sto_50`: Water potential at 50% stomatal closure [MPa]
- `Psi_sto_99`: Water potential at 99% stomatal closure [MPa]
- `CT`: Photosynthesis type (1=C3, 2=C4)
- `Vmax`: Maximum carboxylation rate [μmolCO₂/m²/s]
- `DS`: Deactivation entropy [kJ/mol/K]
- `Ha`: Activation energy [kJ/mol]
- `FI`: Quantum yield [-]
- `Oa`: Atmospheric O₂ concentration [μmolO₂/mol]
- `Do`: Empirical coefficient for stomatal response [Pa]
- `a1`: Empirical coefficient for stomatal response [-]
- `go`: Minimum stomatal conductance [mol/m²/s]
- `e_rel`: Age relative efficiency [-]
- `e_relN`: Nitrogen relative efficiency [-]
- `gmes`: Mesophyll conductance [mol/m²/s]
- `rjv`: Ratio of Jmax to Vmax [-]
- `mSl`: Empirical slope parameter [-]
- `Sl`: Entropy parameter [kJ/mol/K]
- `atol`: Absolute tolerance for CO₂ concentration solver

# Returns
- `rs_sun`: Sunlit stomatal resistance [s/m]
- `rs_shd`: Shaded stomatal resistance [s/m]
- `Ci_sun`: Sunlit leaf internal CO₂ concentration [μmolCO₂/mol]
- `Ci_shd`: Shaded leaf internal CO₂ concentration [μmolCO₂/mol]
- `An`: Net assimilation [μmolCO₂/m²/s]
- `Rdark`: Dark respiration [μmolCO₂/m²/s]
- `Lpho`: Light use [W/m²]
- `SIF`: Solar-induced fluorescence [W/m²/sr/μm]
- `DCi`: Change in CO₂ concentration [μmolCO₂/mol]
"""
function canopy_resistance_an_evolution(
    PAR_sun::FT,
    PAR_shd::FT,
    LAI::FT,
    Kopt::FT,
    Knit::FT,
    Fsun::FT,
    Fshd::FT,
    Citm1_sun::FT,
    Citm1_shd::FT,
    Ca::FT,
    ra::FT,
    rb::FT,
    Ts::FT,
    Ta::FT,
    Pre::FT,
    Ds::FT,
    Psi_L::FT,
    Psi_sto_50::FT,
    Psi_sto_99::FT,
    CT::Int,
    Vmax::FT,
    DS::FT,
    Ha::FT,
    FI::FT,
    Oa::FT,
    Do::FT,
    a1::FT,
    go::FT,
    e_rel::FT,
    e_relN::FT,
    gmes::FT,
    rjv::FT,
    mSl::FT,
    Sl::FT,
    atol::FT,
) where {FT<:AbstractFloat}

    # Ensure minimum CO₂ concentration
    Citm1_sun = max(Citm1_sun, FT(200))
    Citm1_shd = max(Citm1_shd, FT(200))

    # Apply age and nitrogen efficiency to Vmax
    Vmax *= e_rel * e_relN

    # Calculate fractions for scaling
    FsunV = (1 - exp(-Kopt * LAI)) / (Kopt * LAI)
    FsunV = FsunV < 0.01 ? zero(FT) : FsunV
    FsunV = min(FsunV, FT(1))
    FshdV = 1 - FsunV

    # Two big leaves with nitrogen extinction
    Can_sun = (1 - exp(-(Kopt + Knit) * LAI)) / (Kopt + Knit)
    Can_shd = (1 - exp(-Knit * LAI)) / Knit - Can_sun

    # Scale maximum Rubisco capacity
    Vmax_sun = FsunV > 0 ? Vmax * Can_sun / (LAI * FsunV) : zero(FT)
    Vmax_shd = Vmax * Can_shd / (LAI * FshdV)

    # Initialize conductance parameters
    go_sun = go
    rb_sun = rb
    go_shd = go
    rb_shd = rb
    gmes_sun = gmes
    gmes_shd = gmes

    # Scale PAR by leaf area fractions
    PAR_sun = Fsun > 0 ? PAR_sun / (LAI * Fsun) : zero(FT)
    PAR_shd = Fshd > 0 ? PAR_shd / (LAI * Fshd) : zero(FT)

    # Process sunlit fraction
    # TODO: refactor by splitting duplicated code in two
    if Fsun > 0
        Ci_sun = find_zero(
            Cc -> co2_concentration(
                Cc,
                PAR_sun,
                Ca,
                ra,
                rb_sun,
                Ts,
                Pre,
                Ds,
                Psi_L,
                Psi_sto_50,
                Psi_sto_99,
                CT,
                Vmax_sun,
                DS,
                Ha,
                FI,
                Oa,
                Do,
                a1,
                go_sun,
                gmes_sun,
                rjv,
            ),
            Citm1_sun;
            atol=atol,
        )

        CiF_sun, An_sun, rc_sun, Rdark_sun, SIF_sun = photosynthesis_biochemical(
            Ci_sun,
            PAR_sun,
            Ca,
            ra,
            rb_sun,
            Ts,
            Pre,
            Ds,
            Psi_L,
            Psi_sto_50,
            Psi_sto_99,
            CT,
            Vmax_sun,
            DS,
            Ha,
            FI,
            Oa,
            Do,
            a1,
            go_sun,
            gmes_sun,
            rjv,
        )
    else
        Ci_sun = zero(FT)
        CiF_sun = zero(FT)
        An_sun = zero(FT)
        Rdark_sun = zero(FT)
        rc_sun = Inf
        SIF_sun = zero(FT)
    end

    # Process shaded fraction
    if Fshd > 0
        Ci_shd = find_zero(
            Cc -> co2_concentration(
                Cc,
                PAR_shd,
                Ca,
                ra,
                rb_shd,
                Ts,
                Pre,
                Ds,
                Psi_L,
                Psi_sto_50,
                Psi_sto_99,
                CT,
                Vmax_shd,
                DS,
                Ha,
                FI,
                Oa,
                Do,
                a1,
                go_shd,
                gmes_shd,
                rjv,
            ),
            Citm1_shd;
            atol=atol,
        )

        CiF_shd, An_shd, rc_shd, Rdark_shd, SIF_shd = photosynthesis_biochemical(
            Ci_shd,
            PAR_shd,
            Ca,
            ra,
            rb_shd,
            Ts,
            Pre,
            Ds,
            Psi_L,
            Psi_sto_50,
            Psi_sto_99,
            CT,
            Vmax_shd,
            DS,
            Ha,
            FI,
            Oa,
            Do,
            a1,
            go_shd,
            gmes_shd,
            rjv,
        )
    else
        Ci_shd = zero(FT)
        CiF_shd = zero(FT)
        An_shd = zero(FT)
        Rdark_shd = zero(FT)
        rc_shd = Inf
        SIF_shd = zero(FT)
    end

    # Calculate CO₂ concentration changes
    DCi_sun = Ci_sun - CiF_sun
    DCi_shd = Ci_shd - CiF_shd
    DCi = (DCi_sun + DCi_shd) / 2

    # Scale to canopy level
    An = An_sun * (LAI * Fsun) + An_shd * (LAI * Fshd)
    Rdark = Rdark_sun * (LAI * Fsun) + Rdark_shd * (LAI * Fshd)
    SIF = SIF_sun * (LAI * Fsun) + SIF_shd * (LAI * Fshd)

    rs_sun = rc_sun
    rs_shd = rc_shd

    # Calculate light use
    lanp = FT(0.469)  # [J/μmol CO₂] from Blanken et al 1997
    Lpho = (An + Rdark) * lanp  # [W/m²]

    return rs_sun, rs_shd, Ci_sun, Ci_shd, An, Rdark, Lpho, SIF, DCi
end
