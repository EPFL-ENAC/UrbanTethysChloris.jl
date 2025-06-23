"""
    swr_abs_indoors(
        SWRinWsun::FT,
        SWRinWshd::FT,
        Hbuild::FT,
        Wroof::FT,
        abc::FT,
        abw::FT,
        abg::FT,
        abm::FT
    ) where {FT<:AbstractFloat}

Calculate shortwave radiation absorption inside building.

# Arguments
- `SWRinWsun`: Incoming shortwave radiation on sunlit wall [W/m²]
- `SWRinWshd`: Incoming shortwave radiation on shaded wall [W/m²]
- `Hbuild`: Building height [m]
- `Wroof`: Roof width [m]
- `abc`: Ceiling albedo
- `abw`: Wall albedo
- `abg`: Ground albedo
- `abm`: Internal mass albedo

# Returns
- `SWRinB`: Incoming shortwave radiation for each surface [W/m²]
- `SWRoutB`: Outgoing shortwave radiation for each surface [W/m²]
- `SWRabsB`: Absorbed shortwave radiation for each surface [W/m²]
"""
function swr_abs_indoors(
    SWRinWsun::FT, SWRinWshd::FT, Hbuild::FT, Wroof::FT, abc::FT, abw::FT, abg::FT, abm::FT
) where {FT<:AbstractFloat}
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

    # Calculate SWR for sunlit and shaded halves
    SWRinB_wsun, SWRoutB_wsun, SWRabsB_wsun = swr_abs_building_half(
        A_c,
        A_g,
        A_h,
        SWRinWsun,
        F_gc,
        F_gw,
        F_ww,
        F_wg,
        F_wc,
        F_cg,
        F_cw,
        abc,
        abw,
        abm,
        abg,
    )

    SWRinB_wshd, SWRoutB_wshd, SWRabsB_wshd = swr_abs_building_half(
        A_c,
        A_g,
        A_h,
        SWRinWshd,
        F_gc,
        F_gw,
        F_ww,
        F_wg,
        F_wc,
        F_cg,
        F_cw,
        abc,
        abw,
        abm,
        abg,
    )

    # Combine results for whole building
    SWRinB = (
        SWRinCeiling=(SWRinB_wsun.SWRinCeiling + SWRinB_wshd.SWRinCeiling) / 2,
        SWRinWallsun=SWRinB_wsun.SWRinWall,
        SWRinWallshd=SWRinB_wshd.SWRinWall,
        SWRinInternalMass=(SWRinB_wsun.SWRinInternalMass + SWRinB_wsun.SWRinInternalMass),
        SWRinGround=(SWRinB_wsun.SWRinGround + SWRinB_wshd.SWRinGround) / 2,
    )

    SWRoutB = (
        SWRoutCeiling=(SWRoutB_wsun.SWRoutCeiling + SWRoutB_wshd.SWRoutCeiling) / 2,
        SWRoutWallsun=SWRoutB_wsun.SWRoutWall,
        SWRoutWallshd=SWRoutB_wshd.SWRoutWall,
        SWRoutInternalMass=(
            SWRoutB_wsun.SWRoutInternalMass + SWRoutB_wshd.SWRoutInternalMass
        ),
        SWRoutGround=(SWRoutB_wsun.SWRoutGround + SWRoutB_wshd.SWRoutGround) / 2,
    )

    SWRabsB = (
        SWRabsCeiling=(SWRabsB_wsun.SWRabsCeiling + SWRabsB_wshd.SWRabsCeiling) / 2,
        SWRabsWallsun=SWRabsB_wsun.SWRabsWall,
        SWRabsWallshd=SWRabsB_wshd.SWRabsWall,
        SWRabsInternalMass=(
            SWRabsB_wsun.SWRabsInternalMass + SWRabsB_wshd.SWRabsInternalMass
        ),
        SWRabsGround=(SWRabsB_wsun.SWRabsGround + SWRabsB_wshd.SWRabsGround) / 2,
    )

    # Energy balance check
    SWREBinternal =
        Hbuild * (SWRinWsun + SWRinWshd) - Wroof * SWRabsB.SWRabsCeiling -
        Hbuild *
        (SWRabsB.SWRabsWallsun + SWRabsB.SWRabsWallshd + SWRabsB.SWRabsInternalMass) -
        Wroof * SWRabsB.SWRabsGround

    if abs(SWREBinternal) >= FT(1e-6)
        @warn "Building interior shortwave radiation balance is not 0. Please check swr_abs_indoors.jl"
    end

    return SWRinB, SWRoutB, SWRabsB
end
