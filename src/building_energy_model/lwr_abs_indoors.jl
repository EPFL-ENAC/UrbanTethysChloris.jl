"""
    lwr_abs_indoors(
        Tinwallsun::FT,
        Tinwallshd::FT,
        Tceiling::FT,
        Tground::FT,
        Tintmass::FT,
        Hbuild::FT,
        Wroof::FT,
        ec::FT,
        eg::FT,
        em::FT,
        ew::FT
    ) where {FT<:AbstractFloat}

Calculate longwave radiation absorption inside building considering both sunlit and shaded walls.

# Arguments
- `Tinwallsun`: Temperature of sunlit wall [K]
- `Tinwallshd`: Temperature of shaded wall [K]
- `Tceiling`: Ceiling temperature [K]
- `Tground`: Ground temperature [K]
- `Tintmass`: Internal mass temperature [K]
- `Hbuild`: Building height [m]
- `Wroof`: Roof width [m]
- `ec`: Ceiling emissivity
- `eg`: Ground emissivity
- `em`: Internal mass emissivity
- `ew`: Wall emissivity

# Returns
- `LWRinB`: Incoming longwave radiation for each surface [W/m²]
- `LWRoutB`: Outgoing longwave radiation for each surface [W/m²]
- `LWRabsB`: Absorbed longwave radiation for each surface [W/m²]
"""
function lwr_abs_indoors(
    Tinwallsun::FT,
    Tinwallshd::FT,
    Tceiling::FT,
    Tground::FT,
    Tintmass::FT,
    Hbuild::FT,
    Wroof::FT,
    ec::FT,
    eg::FT,
    em::FT,
    ew::FT,
) where {FT<:AbstractFloat}
    @warn "Using untested lwr_abs_building_half function."
    # Half roof width for internal calculations
    Wroofint = Wroof / 2

    # Calculate view factors
    F_gc, F_gw, F_ww, F_wg, F_wc, F_cg, F_cw, _ = view_factor_internal(Hbuild, Wroofint)

    # Normalized surface areas
    A_c = one(FT)
    A_g = one(FT)
    A_h = Hbuild / Wroofint

    # Verify view factors sum to 1
    SVF = [F_gc + 2*F_gw, F_ww + F_wg + F_wc, F_cg + 2*F_cw]
    SVF2 = [F_gc + 2*F_wc*A_h, F_cg + 2*F_wg*A_h, F_ww + F_wc/A_h + F_gw/A_h]

    for s in SVF
        if !(0.999 ≤ s ≤ 1.001)
            @warn "The view factors do not add up to 1 for a canyon with trees"
        end
    end

    # Calculate LWR for sunlit and shaded halves
    LWRinB_wsun, LWRoutB_wsun, LWRabsB_wsun = lwr_abs_building_half(
        Tceiling,
        Tinwallsun,
        Tintmass,
        Tground,
        A_c,
        A_g,
        A_h,
        F_gc,
        F_gw,
        F_ww,
        F_wg,
        F_wc,
        F_cg,
        F_cw,
        ec,
        eg,
        ew,
        em,
    )

    LWRinB_wshd, LWRoutB_wshd, LWRabsB_wshd = lwr_abs_building_half(
        Tceiling,
        Tinwallshd,
        Tintmass,
        Tground,
        A_c,
        A_g,
        A_h,
        F_gc,
        F_gw,
        F_ww,
        F_wg,
        F_wc,
        F_cg,
        F_cw,
        ec,
        eg,
        ew,
        em,
    )

    # Combine results for whole building
    LWRinB = (
        LWRinCeiling=(LWRinB_wsun.LWRinCeiling + LWRinB_wshd.LWRinCeiling) / 2,
        LWRinWallsun=LWRinB_wsun.LWRinWall,
        LWRinWallshd=LWRinB_wshd.LWRinWall,
        LWRinInternalMass=(LWRinB_wsun.LWRinInternalMass + LWRinB_wshd.LWRinInternalMass),
        LWRinGround=(LWRinB_wsun.LWRinGround + LWRinB_wshd.LWRinGround) / 2,
    )

    LWRoutB = (
        LWRoutCeiling=(LWRoutB_wsun.LWRoutCeiling + LWRoutB_wshd.LWRoutCeiling) / 2,
        LWRoutWallsun=LWRoutB_wsun.LWRoutWall,
        LWRoutWallshd=LWRoutB_wshd.LWRoutWall,
        LWRoutInternalMass=(
            LWRoutB_wsun.LWRoutInternalMass + LWRoutB_wshd.LWRoutInternalMass
        ),
        LWRoutGround=(LWRoutB_wsun.LWRoutGround + LWRoutB_wshd.LWRoutGround) / 2,
    )

    LWRabsB = (
        LWRabsCeiling=(LWRabsB_wsun.LWRabsCeiling + LWRabsB_wshd.LWRabsCeiling) / 2,
        LWRabsWallsun=LWRabsB_wsun.LWRabsWall,
        LWRabsWallshd=LWRabsB_wshd.LWRabsWall,
        LWRabsInternalMass=(
            LWRabsB_wsun.LWRabsInternalMass + LWRabsB_wshd.LWRabsInternalMass
        ),
        LWRabsGround=(LWRabsB_wsun.LWRabsGround + LWRabsB_wshd.LWRabsGround) / 2,
    )

    return LWRinB, LWRoutB, LWRabsB
end
