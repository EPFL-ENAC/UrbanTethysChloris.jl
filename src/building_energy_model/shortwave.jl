"""
    swr_abs_building_half(
        A_c::FT,
        A_g::FT,
        A_h::FT,
        SWRinW::FT,
        F_gc::FT,
        F_gw::FT,
        F_ww::FT,
        F_wg::FT,
        F_wc::FT,
        F_cg::FT,
        F_cw::FT,
        abc::FT,
        abw::FT,
        abm::FT,
        abg::FT
    ) where {FT<:AbstractFloat}

Calculate shortwave radiation absorption in building for half of the building.

# Arguments
- `A_c`: Ceiling area [m²]
- `A_g`: Ground area [m²]
- `A_h`: Wall area [m²]
- `SWRinW`: Incoming shortwave radiation on wall [W/m²]
- `F_gc`: View factor ground to ceiling
- `F_gw`: View factor ground to wall
- `F_ww`: View factor wall to wall
- `F_wg`: View factor wall to ground
- `F_wc`: View factor wall to ceiling
- `F_cg`: View factor ceiling to ground
- `F_cw`: View factor ceiling to wall
- `abc`: Ceiling albedo
- `abw`: Wall albedo
- `abm`: Internal mass albedo
- `abg`: Ground albedo

# Returns
- `SWRinB::NamedTuple`: Incoming shortwave radiation for each surface [W/m²]
- `SWRoutB::NamedTuple`: Outgoing shortwave radiation for each surface [W/m²]
- `SWRabsB::NamedTuple`: Absorbed shortwave radiation for each surface [W/m²]
"""
function swr_abs_building_half(
    A_c::FT,
    A_g::FT,
    A_h::FT,
    SWRinW::FT,
    F_gc::FT,
    F_gw::FT,
    F_ww::FT,
    F_wg::FT,
    F_wc::FT,
    F_cg::FT,
    F_cw::FT,
    abc::FT,
    abw::FT,
    abm::FT,
    abg::FT,
) where {FT<:AbstractFloat}
    # Set up albedo vector
    ai = [abc, abw, abm, abg]

    # View factor matrix for infinite reflections
    Tij = [
        1 -abc*F_cw -abc*F_cw -abc*F_cg
        -abw*F_wc 1 -abw*F_ww -abw*F_wg
        -abm*F_wc -abm*F_ww 1 -abm*F_wc
        -abg*F_gc -abg*F_gw -abg*F_gw 1
    ]

    # Incoming shortwave radiation from sky
    Omega_i = [
        abc * F_cw * SWRinW
        zero(FT)
        abm * F_ww * SWRinW
        abg * F_gw * SWRinW
    ]

    # Calculate outgoing radiation
    B_i = Tij \ Omega_i

    # View factor matrix for incoming radiation
    Tij2 = [
        0 F_cw F_cw F_cg
        F_wc 0 F_ww F_wg
        F_wc F_ww 0 F_wc
        F_gc F_gw F_gw 0
    ]

    # Direct radiation
    SWRdir_i = [
        F_cw * SWRinW
        zero(FT)
        F_ww * SWRinW
        F_gw * SWRinW
    ]

    # Calculate incoming radiation
    A_i1 = Tij2 * B_i + SWRdir_i
    A_i = similar(B_i)

    # Handle zero albedo cases
    for i in eachindex(ai)
        if ai[i] ≈ zero(FT)
            A_i[i] = A_i1[i]
        else
            A_i[i] = B_i[i] / ai[i]
        end
    end

    # Calculate net absorbed radiation
    Qnet_i = A_i - B_i

    # Store results
    SWRout_i = B_i
    SWRin_i = A_i
    SWRnet_i = Qnet_i

    # Energy balance checks
    TotalSWRSurface_in =
        (SWRin_i[1]*A_c + SWRin_i[2]*A_h + SWRin_i[3]*A_h + SWRin_i[4]*A_g) / A_g
    TotalSWRSurface_abs =
        (SWRnet_i[1]*A_c + SWRnet_i[2]*A_h + SWRnet_i[3]*A_h + SWRnet_i[4]*A_g) / A_g
    TotalSWRSurface_out =
        (SWRout_i[1]*A_c + SWRout_i[2]*A_h + SWRout_i[3]*A_h + SWRout_i[4]*A_g) / A_g

    EBSurface = TotalSWRSurface_in - TotalSWRSurface_abs - TotalSWRSurface_out
    EBindoors = A_h/A_g * SWRinW - TotalSWRSurface_abs

    if abs(EBindoors) >= FT(1e-6)
        @warn "EBindoors is not 0. Please check swr_abs_building_half.jl"
    end
    if abs(EBSurface) >= FT(1e-6)
        @warn "EBSurface is not 0. Please check swr_abs_building_half.jl"
    end

    # Create return structs
    SWRinB = (
        SWRinCeiling=SWRin_i[1],
        SWRinWall=SWRin_i[2],
        SWRinInternalMass=SWRin_i[3],
        SWRinGround=SWRin_i[4],
    )

    SWRoutB = (
        SWRoutCeiling=SWRout_i[1],
        SWRoutWall=SWRout_i[2],
        SWRoutInternalMass=SWRout_i[3],
        SWRoutGround=SWRout_i[4],
    )

    SWRabsB = (
        SWRabsCeiling=SWRnet_i[1],
        SWRabsWall=SWRnet_i[2],
        SWRabsInternalMass=SWRnet_i[3],
        SWRabsGround=SWRnet_i[4],
    )

    return SWRinB, SWRoutB, SWRabsB
end

"""
    swr_abs_indoors_no_int_mass(
        SWRinWsun::FT,
        SWRinWshd::FT,
        Hbuild::FT,
        Wroof::FT,
        abc::FT,
        abw::FT,
        abg::FT
    ) where {FT<:AbstractFloat}

Calculate shortwave radiation absorption inside building without internal mass.

# Arguments
- `SWRinWsun`: Incoming shortwave radiation on sunlit wall [W/m²]
- `SWRinWshd`: Incoming shortwave radiation on shaded wall [W/m²]
- `Hbuild`: Building height [m]
- `Wroof`: Roof width [m]
- `abc`: Ceiling albedo [-]
- `abw`: Wall albedo [-]
- `abg`: Ground albedo [-]

# Returns
- `SWRinB::NamedTuple`: Incoming shortwave radiation for each surface [W/m²]
- `SWRoutB::NamedTuple`: Outgoing shortwave radiation for each surface [W/m²]
- `SWRabsB::NamedTuple`: Absorbed shortwave radiation for each surface [W/m²]
- `SWREBB::NamedTuple`: Energy balance for each surface [W/m²]
"""
function swr_abs_indoors_no_int_mass(
    SWRinWsun::FT, SWRinWshd::FT, Hbuild::FT, Wroof::FT, abc::FT, abw::FT, abg::FT
) where {FT<:AbstractFloat}
    # Calculate view factors
    F_gc, F_gw, F_ww, F_wg, F_wc, F_cg, F_cw, _ = view_factor_internal(Hbuild, Wroof)

    # Normalized surface areas
    A_c = one(FT)
    A_g = one(FT)
    A_h = Hbuild / Wroof

    # Verify view factors sum to 1
    SVF = [F_gc + 2*F_gw, F_ww + F_wg + F_wc, F_cg + 2*F_cw]
    SVF2 = [F_gc + 2*F_wc*A_h, F_cg + 2*F_wg*A_h, F_ww + F_wc/A_h + F_gw/A_h]

    for s in SVF
        if !(0.999 ≤ s ≤ 1.001)
            @warn "The view factors do not add up to 1 for a canyon with trees"
        end
    end

    # Set up albedo vector
    ai = [abc, abw, abw, abg]

    # View factor matrix for infinite reflections
    Tij = [
        1 -abc*F_cw -abc*F_cw -abc*F_cg
        -abw*F_wc 1 -abw*F_ww -abw*F_wg
        -abw*F_wc -abw*F_ww 1 -abw*F_wc
        -abg*F_gc -abg*F_gw -abg*F_gw 1
    ]

    # Incoming shortwave radiation from sky
    Omega_i = [
        abc * (F_cw*SWRinWsun + F_cw*SWRinWshd)
        abw * F_ww * SWRinWshd
        abw * F_ww * SWRinWsun
        abg * (F_gw*SWRinWsun + F_gw*SWRinWshd)
    ]

    # Calculate outgoing radiation
    B_i = Tij \ Omega_i

    # Calculate incoming radiation
    Tij2 = [
        0 F_cw F_cw F_cg
        F_wc 0 F_ww F_wg
        F_wc F_ww 0 F_wc
        F_gc F_gw F_gw 0
    ]

    SWRdir_i = [
        F_cw*SWRinWsun + F_cw*SWRinWshd
        F_ww*SWRinWshd
        F_ww*SWRinWsun
        F_gw*SWRinWsun + F_gw*SWRinWshd
    ]

    A_i1 = Tij2 * B_i + SWRdir_i
    A_i = similar(B_i)

    # Calculate incoming radiation accounting for zero albedos
    for i in eachindex(ai)
        if ai[i] ≈ zero(FT)
            A_i[i] = A_i1[i]
        else
            A_i[i] = B_i[i] / ai[i]
        end
    end

    # Calculate net absorbed radiation
    Qnet_i = A_i - B_i

    # Store results
    SWRout_i = B_i
    SWRin_i = A_i
    SWRnet_i = Qnet_i

    # Energy balance checks
    TotalSWRSurface_in =
        (SWRin_i[1]*A_c + SWRin_i[2]*A_h + SWRin_i[3]*A_h + SWRin_i[4]*A_g) / A_g
    TotalSWRSurface_abs =
        (SWRnet_i[1]*A_c + SWRnet_i[2]*A_h + SWRnet_i[3]*A_h + SWRnet_i[4]*A_g) / A_g
    TotalSWRSurface_out =
        (SWRout_i[1]*A_c + SWRout_i[2]*A_h + SWRout_i[3]*A_h + SWRout_i[4]*A_g) / A_g

    EBSurface = TotalSWRSurface_in - TotalSWRSurface_abs - TotalSWRSurface_out
    EBindoors = A_h/A_g * (SWRinWsun + SWRinWshd) - TotalSWRSurface_abs

    if abs(EBindoors) >= FT(1e-6)
        @warn "EBindoors is not 0. Please check swr_abs_indoors_no_int_mass.jl"
    end
    if abs(EBSurface) >= FT(1e-6)
        @warn "EBSurface is not 0. Please check swr_abs_indoors_no_int_mass.jl"
    end

    # Create return structs
    SWRinB = (
        SWRinCeiling=SWRin_i[1],
        SWRinWallsun=SWRin_i[2],
        SWRinWallshd=SWRin_i[3],
        SWRinGround=SWRin_i[4],
    )

    SWRoutB = (
        SWRoutCeiling=SWRout_i[1],
        SWRoutWallsun=SWRout_i[2],
        SWRoutWallshd=SWRout_i[3],
        SWRoutGround=SWRout_i[4],
    )

    SWRabsB = (
        SWRabsCeiling=SWRnet_i[1],
        SWRabsWallsun=SWRnet_i[2],
        SWRabsWallshd=SWRnet_i[3],
        SWRabsGround=SWRnet_i[4],
    )

    SWREBB = (
        SWREBCeiling=SWRinB.SWRinCeiling - SWRoutB.SWRoutCeiling - SWRabsB.SWRabsCeiling,
        SWREBWallsun=SWRinB.SWRinWallsun - SWRoutB.SWRoutWallsun - SWRabsB.SWRabsWallsun,
        SWREBWallshd=SWRinB.SWRinWallshd - SWRoutB.SWRoutWallshd - SWRabsB.SWRabsWallshd,
        SWREBGround=SWRinB.SWRinGround - SWRoutB.SWRoutGround - SWRabsB.SWRabsGround,
    )

    # Check energy balance for each surface
    for (component, val) in pairs(SWREBB)
        if abs(val) >= FT(1e-6)
            @warn "$(component) is not 0. Please check swr_abs_indoors_no_int_mass.jl"
        end
    end

    return SWRinB, SWRoutB, SWRabsB, SWREBB
end

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
- `SWRinB::NamedTuple`: Incoming shortwave radiation for each surface [W/m²]
- `SWRoutB::NamedTuple`: Outgoing shortwave radiation for each surface [W/m²]
- `SWRabsB::NamedTuple`: Absorbed shortwave radiation for each surface [W/m²]
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
