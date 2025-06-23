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
- `SWRinB`: Incoming shortwave radiation for each surface [W/m²]
- `SWRoutB`: Outgoing shortwave radiation for each surface [W/m²]
- `SWRabsB`: Absorbed shortwave radiation for each surface [W/m²]
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
