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
- `LWRinB`: Incoming longwave radiation for each surface [W/m²]
- `LWRoutB`: Outgoing longwave radiation for each surface [W/m²]
- `LWRabsB`: Absorbed longwave radiation for each surface [W/m²]
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

    # Energy balance check
    TotalLWRSurface_in =
        (LWRin_i[1]*A_c + LWRin_i[2]*A_h + LWRin_i[3]*A_h + LWRin_i[4]*A_g) / A_g
    TotalLWRSurface_abs =
        (Qnet_i[1]*A_c + Qnet_i[2]*A_h + Qnet_i[3]*A_h + Qnet_i[4]*A_g) / A_g
    TotalLWRSurface_out = (B_i[1]*A_c + B_i[2]*A_h + B_i[3]*A_h + B_i[4]*A_g) / A_g

    EBSurface = TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out
    if abs(EBSurface) >= FT(1e-6)
        @warn "EBSurface is not 0. Please check lwr_abs_building_half.jl"
    end

    # Create return structs
    LWRinB = (
        LWRinCeiling=A_i[1], LWRinWall=A_i[2], LWRinInternalMass=A_i[3], LWRinGround=A_i[4]
    )

    LWRoutB = (
        LWRoutCeiling=B_i[1],
        LWRoutWall=B_i[2],
        LWRoutInternalMass=B_i[3],
        LWRoutGround=B_i[4],
    )

    LWRabsB = (
        LWRabsCeiling=Qnet_i[1],
        LWRabsWall=Qnet_i[2],
        LWRabsInternalMass=Qnet_i[3],
        LWRabsGround=Qnet_i[4],
    )

    return LWRinB, LWRoutB, LWRabsB
end
