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
- `LWRinB`: Incoming longwave radiation for each surface [W/m²]
- `LWRoutB`: Outgoing longwave radiation for each surface [W/m²]
- `LWRabsB`: Absorbed longwave radiation for each surface [W/m²]
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

    return LWRinB, LWRoutB, LWRabsB
end
