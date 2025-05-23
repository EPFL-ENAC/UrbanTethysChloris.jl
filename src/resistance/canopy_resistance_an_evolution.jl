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
    Ci_sun, CiF_sun, An_sun, rc_sun, Rdark_sun, SIF_sun = process_leaf_fraction(
        Fsun,
        PAR_sun,
        Citm1_sun,
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
        atol,
    )

    # Process shaded fraction
    Ci_shd, CiF_shd, An_shd, rc_shd, Rdark_shd, SIF_shd = process_leaf_fraction(
        Fshd,
        PAR_shd,
        Citm1_shd,
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
        atol,
    )

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

"""
    process_leaf_fraction(
        F::FT,
        PAR::FT,
        Citm1::FT,
        Ca::FT,
        ra::FT,
        rb::FT,
        Ts::FT,
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
        gmes::FT,
        rjv::FT,
        atol::FT
    ) where {FT<:AbstractFloat}

Process CO₂ assimilation for a leaf fraction (sunlit or shaded).

# Arguments
- `F`: Leaf fraction [-]
- `PAR`: Photosynthetically active radiation [W/m²]
- `Citm1`: Previous CO₂ concentration [μmolCO₂/mol]
- `Ca`: Atmospheric CO₂ concentration [μmolCO₂/mol]
- `ra`: Aerodynamic resistance [s/m]
- `rb`: Boundary layer resistance [s/m]
- `Ts`: Leaf temperature [K]
- `Pre`: Atmospheric pressure [Pa]
- `Ds`: Vapor pressure deficit [Pa]
- `Psi_L`: Leaf water potential [MPa]
- `Psi_sto_50`: 50% storage water potential [MPa]
- `Psi_sto_99`: 99% storage water potential [MPa]
- `CT`: Temperature coefficient
- `Vmax`: Maximum carboxylation rate [μmolCO₂/m²/s]
- `DS`: Stomatal conductance [m/s]
- `Ha`: Leaf area [m²]
- `FI`: Light intensity [μmol photons/m²/s]
- `Oa`: Ambient oxygen concentration [μmolO₂/mol]
- `Do`: Diffusion coefficient [m²/s]
- `a1`: Parameter for photosynthesis model
- `go`: Stomatal conductance [m/s]
- `gmes`: Mesophyll conductance [m/s]
- `rjv`: Jmax/Rubisco parameter
- `atol`: Absolute tolerance for root finding

# Returns
- `Ci`: Leaf internal CO₂ concentration [μmolCO₂/mol]
- `CiF`: Leaf internal CO₂ concentration at equilibrium [μmolCO₂/mol]
- `An`: Net assimilation rate [μmolCO₂/m²/s]
- `rc`: Canopy resistance [s/m]
- `Rdark`: Dark respiration rate [μmolCO₂/m²/s]
- `SIF`: Solar-induced fluorescence [W/m²/sr/μm]
"""
function process_leaf_fraction(
    F::FT,
    PAR::FT,
    Citm1::FT,
    Ca::FT,
    ra::FT,
    rb::FT,
    Ts::FT,
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
    gmes::FT,
    rjv::FT,
    atol::FT,
) where {FT<:AbstractFloat}
    if F > 0
        Ci = find_zero(
            Cc -> co2_concentration(
                Cc,
                PAR,
                Ca,
                ra,
                rb,
                Ts,
                Pre,
                Ds,
                Psi_L,
                Psi_sto_50,
                Psi_sto_99,
                CT,
                Vmax,
                DS,
                Ha,
                FI,
                Oa,
                Do,
                a1,
                go,
                gmes,
                rjv,
            ),
            Citm1;
            atol=atol,
        )

        CiF, An, rc, Rdark, SIF = photosynthesis_biochemical(
            Ci,
            PAR,
            Ca,
            ra,
            rb,
            Ts,
            Pre,
            Ds,
            Psi_L,
            Psi_sto_50,
            Psi_sto_99,
            CT,
            Vmax,
            DS,
            Ha,
            FI,
            Oa,
            Do,
            a1,
            go,
            gmes,
            rjv,
        )
    else
        Ci = zero(FT)
        CiF = zero(FT)
        An = zero(FT)
        rc = Inf
        Rdark = zero(FT)
        SIF = zero(FT)
    end

    return Ci, CiF, An, rc, Rdark, SIF
end
