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
- `SWRinB::FT`: Incoming shortwave radiation for each surface [W/m²]
- `SWRoutB::FT`: Outgoing shortwave radiation for each surface [W/m²]
- `SWRabsB::FT`: Absorbed shortwave radiation for each surface [W/m²]
- `SWREBB::FT`: Energy balance for each surface [W/m²]
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
