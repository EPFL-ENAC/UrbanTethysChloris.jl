"""
    lwr_abs_building_half(
        Tceiling::FT,
        Tinwall::FT,
        Tintmass::FT,
        Tground::FT,
        A_c::FT,
        A_g::FT,
        A_h::FT,
        F_gc::FT,
        F_gw::FT,
        F_ww::FT,
        F_wg::FT,
        F_wc::FT,
        F_cg::FT,
        F_cw::FT,
        ec::FT,
        eg::FT,
        ew::FT,
        em::FT
    ) where {FT<:AbstractFloat}

Calculate longwave radiation absorption in building for half of the building.

# Arguments
- `Tceiling`: Ceiling temperature [K]
- `Tinwall`: Interior wall temperature [K]
- `Tintmass`: Internal mass temperature [K]
- `Tground`: Ground temperature [K]
- `A_c`: Ceiling area [m²]
- `A_g`: Ground area [m²]
- `A_h`: Wall area [m²]
- `F_gc`: View factor ground to ceiling
- `F_gw`: View factor ground to wall
- `F_ww`: View factor wall to wall
- `F_wg`: View factor wall to ground
- `F_wc`: View factor wall to ceiling
- `F_cg`: View factor ceiling to ground
- `F_cw`: View factor ceiling to wall
- `ec`: Ceiling emissivity
- `eg`: Ground emissivity
- `ew`: Wall emissivity
- `em`: Internal mass emissivity

# Returns
- `LWRinB::NamedTuple`: Incoming longwave radiation for each surface [W/m²]
- `LWRoutB::NamedTuple`: Outgoing longwave radiation for each surface [W/m²]
- `LWRabsB::NamedTuple`: Absorbed longwave radiation for each surface [W/m²]
"""
function lwr_abs_building_half(
    Tceiling::FT,
    Tinwall::FT,
    Tintmass::FT,
    Tground::FT,
    A_c::FT,
    A_g::FT,
    A_h::FT,
    F_gc::FT,
    F_gw::FT,
    F_ww::FT,
    F_wg::FT,
    F_wc::FT,
    F_cg::FT,
    F_cw::FT,
    ec::FT,
    eg::FT,
    ew::FT,
    em::FT,
) where {FT<:AbstractFloat}
    @warn "Using untested lwr_abs_building_half function."
    # Constants
    bolzm = FT(5.67e-8)  # Stefan-Boltzmann constant [W*m^-2*K-4]

    # View factor matrix for infinite reflections
    Tij = [
        1 -(1-ec)*F_cw -(1-ec)*F_cw -(1-ec)*F_cg
        -(1-ew)*F_wc 1 -(1-ew)*F_ww -(1-ew)*F_wg
        -(1-em)*F_wc -(1-em)*F_ww 1 -(1-em)*F_wg
        -(1-eg)*F_gc -(1-eg)*F_gw -(1-eg)*F_gw 1
    ]

    # Emitted radiation per surface
    Omega_i = [
        ec * bolzm * Tceiling^4
        ew * bolzm * Tinwall^4
        em * bolzm * Tintmass^4
        eg * bolzm * Tground^4
    ]

    # Calculate outgoing radiation
    B_i = Tij \ Omega_i

    # View factor matrix for incoming radiation
    Tij2 = [
        0 F_cw F_cw F_cg
        F_wc 0 F_ww F_wg
        F_wc F_ww 0 F_wg
        F_gc F_gw F_gw 0
    ]

    # Calculate incoming and net radiation
    A_i = Tij2 * B_i
    e_i = [ec, ew, em, eg]

    # Calculate net absorbed radiation
    Qnet_i = similar(B_i)
    for i in eachindex(e_i)
        if e_i[i] == 1
            Qnet_i[i] = A_i[i] - Omega_i[i]
        else
            Qnet_i[i] = (e_i[i] * B_i[i] - Omega_i[i]) / (1 - e_i[i])
        end
    end

    LWRout_i = B_i
    LWRemit_i = Omega_i
    LWRin_i = A_i
    LWRnet_i = Qnet_i

    # Energy balance check
    TotalLWRSurface_in =
        (LWRin_i[1]*A_c + LWRin_i[2]*A_h + LWRin_i[3]*A_h + LWRin_i[4]*A_g) / A_g
    TotalLWRSurface_abs =
        (LWRnet_i[1]*A_c + LWRnet_i[2]*A_h + LWRnet_i[3]*A_h + LWRnet_i[4]*A_g) / A_g
    TotalLWRSurface_out =
        (LWRout_i[1]*A_c + LWRout_i[2]*A_h + LWRout_i[3]*A_h + LWRout_i[4]*A_g) / A_g

    EBSurface = TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out
    if abs(EBSurface) >= FT(1e-6)
        @warn "EBSurface is not 0. Please check lwr_abs_building_half.jl"
    end

    # Create return structs
    LWRinB = (
        LWRinCeiling=LWRin_i[1],
        LWRinWall=LWRin_i[2],
        LWRinInternalMass=LWRin_i[3],
        LWRinGround=LWRin_i[4],
    )

    LWRoutB = (
        LWRoutCeiling=LWRout_i[1],
        LWRoutWall=LWRout_i[2],
        LWRoutInternalMass=LWRout_i[3],
        LWRoutGround=LWRout_i[4],
    )

    LWRabsB = (
        LWRabsCeiling=LWRnet_i[1],
        LWRabsWall=LWRnet_i[2],
        LWRabsInternalMass=LWRnet_i[3],
        LWRabsGround=LWRnet_i[4],
    )

    return LWRinB, LWRoutB, LWRabsB
end

"""
    lwr_abs_indoors_no_int_mass(
        Tinwallsun::FT,
        Tinwallshd::FT,
        Tceiling::FT,
        Tground::FT,
        Hbuild::FT,
        Wroof::FT,
        ec::FT,
        eg::FT,
        ew::FT
    ) where {FT<:AbstractFloat}

Calculate longwave radiation absorption inside building without internal mass.

# Arguments
- `Tinwallsun`: Temperature of sunlit wall [K]
- `Tinwallshd`: Temperature of shaded wall [K]
- `Tceiling`: Ceiling temperature [K]
- `Tground`: Ground temperature [K]
- `Hbuild`: Building height [m]
- `Wroof`: Roof width [m]
- `ec`: Ceiling emissivity
- `eg`: Ground emissivity
- `ew`: Wall emissivity

# Returns
- `LWRinB::NamedTuple`: Incoming longwave radiation for each surface [W/m²]
- `LWRoutB::NamedTuple`: Outgoing longwave radiation for each surface [W/m²]
- `LWRabsB::NamedTuple`: Absorbed longwave radiation for each surface [W/m²]
"""
function lwr_abs_indoors_no_int_mass(
    Tinwallsun::FT,
    Tinwallshd::FT,
    Tceiling::FT,
    Tground::FT,
    Hbuild::FT,
    Wroof::FT,
    ec::FT,
    eg::FT,
    ew::FT,
) where {FT<:AbstractFloat}
    # Calculate view factors
    F_gc, F_gw, F_ww, F_wg, F_wc, F_cg, F_cw, _ = view_factor_internal(Hbuild, Wroof)

    # Normalized surface areas
    A_c = Wroof / Wroof
    A_g = Wroof / Wroof
    A_h = Hbuild / Wroof
    bolzm = FT(5.67e-8)  # Stefan-Boltzmann constant [W/m²K⁴]

    # Verify view factors sum to 1
    SVF = [F_gc + 2*F_gw, F_ww + F_wg + F_wc, F_cg + 2*F_cw]
    SVF2 = [F_gc + 2*F_wc*A_h, F_cg + 2*F_wg*A_h, F_ww + F_wc/A_h + F_gw/A_h]

    for s in SVF
        if !(0.999 ≤ s ≤ 1.001)
            @warn "The view factors do not add up to 1 for a canyon with trees"
        end
    end
    # TODO: check if call to lwr_abs_building_half instead is possible

    # Calculate infinite reflections
    # View factor matrix for infinite reflections equation
    Tij = [
        1 -(1-ec)*F_cw -(1-ec)*F_cw -(1-ec)*F_cg
        -(1-ew)*F_wc 1 -(1-ew)*F_ww -(1-ew)*F_wg
        -(1-ew)*F_wc -(1-ew)*F_ww 1 -(1-ew)*F_wg
        -(1-eg)*F_gc -(1-eg)*F_gw -(1-eg)*F_gw 1
    ]

    # Emitted radiation per surface
    Omega_i = [
        ec * bolzm * Tceiling^4
        ew * bolzm * Tinwallsun^4
        ew * bolzm * Tinwallshd^4
        eg * bolzm * Tground^4
    ]

    # Calculate outgoing radiation
    B_i = Tij \ Omega_i

    # Calculate incoming radiation
    Tij2 = [
        0 F_cw F_cw F_cg
        F_wc 0 F_ww F_wg
        F_wc F_ww 0 F_wg
        F_gc F_gw F_gw 0
    ]

    # Calculate incoming and net radiation
    A_i = Tij2 * B_i
    e_i = [ec, ew, ew, eg]

    # Calculate net absorbed radiation
    Qnet_i = similar(B_i)
    for i in eachindex(e_i)
        if e_i[i] == 1
            Qnet_i[i] = A_i[i] - Omega_i[i]
        else
            Qnet_i[i] = (e_i[i] * B_i[i] - Omega_i[i]) / (1 - e_i[i])
        end
    end

    # Store results
    LWRout_i = B_i
    LWRin_i = A_i
    LWRnet_i = Qnet_i

    # Energy balance checks
    TotalLWRSurface_in =
        (LWRin_i[1]*A_c + LWRin_i[2]*A_h + LWRin_i[3]*A_h + LWRin_i[4]*A_g) / A_g
    TotalLWRSurface_abs =
        (LWRnet_i[1]*A_c + LWRnet_i[2]*A_h + LWRnet_i[3]*A_h + LWRnet_i[4]*A_g) / A_g
    TotalLWRSurface_out =
        (LWRout_i[1]*A_c + LWRout_i[2]*A_h + LWRout_i[3]*A_h + LWRout_i[4]*A_g) / A_g

    EBSurface = TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out
    if abs(EBSurface) >= FT(1e-6)
        @warn "EBSurface is not 0. Please check lwr_abs_indoors_no_int_mass.jl"
    end

    # Create return structs
    LWRinB = (
        LWRinCeiling=LWRin_i[1],
        LWRinWallsun=LWRin_i[2],
        LWRinWallshd=LWRin_i[3],
        LWRinGround=LWRin_i[4],
    )

    LWRoutB = (
        LWRoutCeiling=LWRout_i[1],
        LWRoutWallsun=LWRout_i[2],
        LWRoutWallshd=LWRout_i[3],
        LWRoutGround=LWRout_i[4],
    )

    LWRabsB = (
        LWRabsCeiling=LWRnet_i[1],
        LWRabsWallsun=LWRnet_i[2],
        LWRabsWallshd=LWRnet_i[3],
        LWRabsGround=LWRnet_i[4],
    )

    # Add LWREBB calculation
    LWREBB = (
        LWREBCeiling=LWRinB.LWRinCeiling - LWRoutB.LWRoutCeiling - LWRabsB.LWRabsCeiling,
        LWREBWallsun=LWRinB.LWRinWallsun - LWRoutB.LWRoutWallsun - LWRabsB.LWRabsWallsun,
        LWREBWallshd=LWRinB.LWRinWallshd - LWRoutB.LWRoutWallshd - LWRabsB.LWRabsWallshd,
        LWREBGround=LWRinB.LWRinGround - LWRoutB.LWRoutGround - LWRabsB.LWRabsGround,
    )

    # Check energy balance for each surface
    for (component, val) in pairs(LWREBB)
        if abs(val) >= FT(1e-6)
            @warn "$(component) is not 0. Please check lwr_abs_indoors_no_int_mass.jl"
        end
    end

    return LWRinB, LWRoutB, LWRabsB, LWREBB
end

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
- `LWRinB::NamedTuple`: Incoming longwave radiation for each surface [W/m²]
- `LWRoutB::NamedTuple`: Outgoing longwave radiation for each surface [W/m²]
- `LWRabsB::NamedTuple`: Absorbed longwave radiation for each surface [W/m²]
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
    @warn "Using untested lwr_abs_indoors function."
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
