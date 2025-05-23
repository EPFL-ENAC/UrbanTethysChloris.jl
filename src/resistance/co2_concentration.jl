"""
    co2_concentration(Cc, IPAR, Csl, ra, rb, Ts, Pre, Ds, Psi_L, Psi_sto_50, Psi_sto_99,
        CT, Vmax, DSE, Ha, FI, Oa, Do, a1, go, gmes, rjv)

Calculate CO₂ concentration difference for root finding in stomatal conductance calculations.

# Arguments
- Same as photosynthesis_biochemical function

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
