"""
    co2_concentration(
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
        Psi_sto_99::FT,
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
        rjv::FT
    ) where {FT<:AbstractFloat}

Calculate CO₂ concentration difference for root finding in stomatal conductance calculations.

# Arguments
- `Cc`: Initial internal CO₂ concentration [μmolCO₂/mol]
- `IPAR`: Absorbed photosynthetically active radiation [W/m²]
- `Csl`: Atmospheric CO₂ concentration [μmolCO₂/mol]
- `ra`: Aerodynamic resistance [s/m]
- `rb`: Leaf boundary layer resistance [s/m]
- `Ts`: Leaf temperature [°C]
- `Pre`: Atmospheric pressure [kPa]
- `Ds`: Leaf-to-air vapor pressure deficit [Pa]
- `Psi_L`: Leaf water potential [MPa]
- `Psi_sto_50`: Leaf water potential at 50% stomatal closure [MPa]
- `Psi_sto_99`: Leaf water potential at 99% stomatal closure [MPa]
- `CT`: Photosynthesis model type (0 = C3, 1 = C4)
- `Vmax`: Maximum carboxylation rate [μmolCO₂/m²s]
- `DSE`: Entropy term [J/mol/K]
- `Ha`: Activation energy [J/mol]
- `FI`: Internal CO₂ to atmospheric CO₂ ratio [-]
- `Oa`: Atmospheric O₂ concentration [mmol/mol]
- `Do`: Empirical coefficient for stomatal sensitivity to VPD [kPa]
- `a1`: Empirical coefficient for stomatal conductance model [-]
- `go`: Minimum stomatal conductance [mol/m²s]
- `gmes`: Mesophyll conductance [mol/m²s]
- `rjv`: Ratio of Jmax to Vmax [-]

# Returns
- `DCi`: Difference between input and calculated CO₂ concentrations [μmolCO₂/mol]
"""
function co2_concentration(
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
    Psi_sto_99::FT,
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
    CcF, = photosynthesis_biochemical(
        Cc,
        IPAR,
        Csl,
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
        DSE,
        Ha,
        FI,
        Oa,
        Do,
        a1,
        go,
        gmes,
        rjv,
    )

    DCi = Cc - CcF

    return DCi
end
