"""
    longwave_absorbed_no_tree(
            h_can::FT, w_can::FT, LWR::FT,
            fgveg::FT, fgbare::FT, fgimp::FT,
            ew::FT, egveg::FT, egbare::FT, egimp::FT,
            Tgimp::FT, Tgbare::FT, Tgveg::FT,
            Twsun::FT, Twshade::FT,
            view_factor::ViewFactor{FT}
    ) where {FT<:AbstractFloat}

Calculate longwave radiation exchange in an urban canyon without trees.

# Arguments
- `h_can`: Building height [-]
- `w_can`: Ground width [-]
- `LWR`: Atmospheric longwave radiation [W/m²]
- `fgveg`: Partitioning ground vegetation [-]
- `fgbare`: Partitioning ground bare [-]
- `fgimp`: Partitioning ground impervious [-]
- `ew`: Wall emissivity [-]
- `egveg`: Ground vegetation emissivity [-]
- `egbare`: Ground bare emissivity [-]
- `egimp`: Ground impervious emissivity [-]
- `Tgimp`: Temperature of ground impervious [K]
- `Tgbare`: Temperature of ground bare [K]
- `Tgveg`: Temperature of ground vegetated [K]
- `Twsun`: Temperature of wall sun [K]
- `Twshade`: Temperature of wall shade [K]
- `view_factor`: View factors for radiation exchange

# Returns
- Tuple of (LWRin_nT, LWRout_nT, LWRabs_nT) containing longwave radiation components
"""
function longwave_absorbed_no_tree(
    h_can::FT,
    w_can::FT,
    LWR::FT,
    fgveg::FT,
    fgbare::FT,
    fgimp::FT,
    ew::FT,
    egveg::FT,
    egbare::FT,
    egimp::FT,
    Tgimp::FT,
    Tgbare::FT,
    Tgveg::FT,
    Twsun::FT,
    Twshade::FT,
    view_factor::RayTracing.ViewFactor{FT},
) where {FT<:AbstractFloat}

    # Stefan-Boltzmann constant [W⋅m⁻²⋅K⁻⁴]
    bolzm = FT(5.67e-8)

    # Check if view factors add up to 1
    SVF = [
        view_factor.F_gs_nT + 2*view_factor.F_gw_nT,
        view_factor.F_ww_nT + view_factor.F_wg_nT + view_factor.F_ws_nT,
        view_factor.F_sg_nT + 2*view_factor.F_sw_nT,
    ]

    all(x -> FT(0.999) ≤ x ≤ FT(1.001), SVF) || @warn "View factors do not add up to 1"

    # Surface indicators
    Cimp = fgimp > 0
    Cbare = fgbare > 0
    Cveg = fgveg > 0

    # View factor matrix for infinite reflections
    Tij = [
        1 0 0 -(1-egveg)*view_factor.F_gw_nT*Cveg -(1-egveg)*view_factor.F_gw_nT*Cveg -(1-egveg)*view_factor.F_gs_nT*Cveg;
        0 1 0 -(1-egbare)*view_factor.F_gw_nT*Cbare -(1-egbare)*view_factor.F_gw_nT*Cbare -(1-egbare)*view_factor.F_gs_nT*Cbare;
        0 0 1 -(1-egimp)*view_factor.F_gw_nT*Cimp -(1-egimp)*view_factor.F_gw_nT*Cimp -(1-egimp)*view_factor.F_gs_nT*Cimp;
        -(1-ew)*view_factor.F_wg_nT*fgveg*Cveg -(1-ew)*view_factor.F_wg_nT*fgbare*Cbare -(1-ew)*view_factor.F_wg_nT*fgimp*Cimp 1 -(1-ew)*view_factor.F_ww_nT -(1-ew)*view_factor.F_ws_nT;
        -(1-ew)*view_factor.F_wg_nT*fgveg*Cveg -(1-ew)*view_factor.F_wg_nT*fgbare*Cbare -(1-ew)*view_factor.F_wg_nT*fgimp*Cimp -(1-ew)*view_factor.F_ww_nT 1 -(1-ew)*view_factor.F_ws_nT;
        0 0 0 0 0 1
    ]

    # Emitted radiation per surface
    Omega_i = [
        egveg*bolzm*(Tgveg^4)*Cveg;
        egbare*bolzm*(Tgbare^4)*Cbare;
        egimp*bolzm*(Tgimp^4)*Cimp;
        ew*bolzm*(Twsun^4);
        ew*bolzm*(Twshade^4);
        LWR
    ]

    # Solve for outgoing radiation
    B_i = Tij \ Omega_i

    # Matrix for incoming radiation
    Tij2 = [
        0 0 0 view_factor.F_gw_nT*Cveg view_factor.F_gw_nT*Cveg view_factor.F_gs_nT*Cveg;
        0 0 0 view_factor.F_gw_nT*Cbare view_factor.F_gw_nT*Cbare view_factor.F_gs_nT*Cbare;
        0 0 0 view_factor.F_gw_nT*Cimp view_factor.F_gw_nT*Cimp view_factor.F_gs_nT*Cimp;
        view_factor.F_wg_nT*fgveg*Cveg view_factor.F_wg_nT*fgbare*Cbare view_factor.F_wg_nT*fgimp*Cimp 0 view_factor.F_ww_nT view_factor.F_ws_nT;
        view_factor.F_wg_nT*fgveg*Cveg view_factor.F_wg_nT*fgbare*Cbare view_factor.F_wg_nT*fgimp*Cimp view_factor.F_ww_nT 0 view_factor.F_ws_nT;
        0 0 0 0 0 0
    ]

    # Calculate incoming and net radiation
    A_i = Tij2 * B_i
    e_i = [egveg; egbare; egimp; ew; ew; FT(0)]

    # Calculate net absorbed radiation
    Qnet_i = similar(B_i)
    for i in eachindex(e_i)
        if e_i[i] == 1
            Qnet_i[i] = A_i[i] - Omega_i[i]
        else
            Qnet_i[i] = (e_i[i]*B_i[i] - Omega_i[i])/(1-e_i[i])
        end
    end
    Qnet_i[6] = 0  # Sky has fixed emission

    # Create output structures
    LWRin_nT = LongwaveRadiation{FT}(;
        GroundImp=A_i[3]*Cimp,
        GroundBare=A_i[2]*Cbare,
        GroundVeg=A_i[1]*Cveg,
        Tree=0,
        WallSun=A_i[4],
        WallShade=A_i[5],
        TotalGround=fgveg*A_i[1] + fgbare*A_i[2] + fgimp*A_i[3],
        TotalCanyon=A_i[1]*fgveg +
                    A_i[2]*fgbare +
                    A_i[3]*fgimp +
                    A_i[4]*h_can/w_can +
                    A_i[5]*h_can/w_can,
    )

    # TODO: create a closure to create the composite structs

    LWRout_nT = LongwaveRadiation{FT}(;
        GroundImp=B_i[3]*Cimp,
        GroundBare=B_i[2]*Cbare,
        GroundVeg=B_i[1]*Cveg,
        Tree=0,
        WallSun=B_i[4],
        WallShade=B_i[5],
        TotalGround=fgveg*B_i[1] + fgbare*B_i[2] + fgimp*B_i[3],
        TotalCanyon=B_i[1]*fgveg +
                    B_i[2]*fgbare +
                    B_i[3]*fgimp +
                    B_i[4]*h_can/w_can +
                    B_i[5]*h_can/w_can,
    )

    LWRabs_nT = LongwaveRadiation{FT}(;
        GroundImp=Qnet_i[3]*Cimp,
        GroundBare=Qnet_i[2]*Cbare,
        GroundVeg=Qnet_i[1]*Cveg,
        Tree=0,
        WallSun=Qnet_i[4],
        WallShade=Qnet_i[5],
        TotalGround=fgveg*Qnet_i[1] + fgbare*Qnet_i[2] + fgimp*Qnet_i[3],
        TotalCanyon=Qnet_i[1]*fgveg +
                    Qnet_i[2]*fgbare +
                    Qnet_i[3]*fgimp +
                    Qnet_i[4]*h_can/w_can +
                    Qnet_i[5]*h_can/w_can,
    )

    LWREB_nT = LWRin_nT - LWRout_nT - LWRabs_nT

    # Energy balance check variables
    A_g = w_can
    A_w = h_can

    TotalLWRSurface_in = LWRin_nT.TotalCanyon
    TotalLWRSurface_abs = LWRabs_nT.TotalCanyon
    TotalLWRSurface_out = LWRout_nT.TotalCanyon
    TotalLWRref_to_atm =
        B_i[1]*view_factor.F_sg_nT*fgveg +
        B_i[2]*view_factor.F_sg_nT*fgbare +
        B_i[3]*view_factor.F_sg_nT*fgimp +
        B_i[4]*view_factor.F_sw_nT +
        B_i[5]*view_factor.F_sw_nT

    # Energy balance checks
    EBSurface = TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out
    EBCanyon = LWR - TotalLWRSurface_abs - TotalLWRref_to_atm

    if abs(EBSurface) >= FT(1e-6)
        @warn "EBSurface is not 0. Please check longwave_absorbed_no_tree"
    end
    if abs(EBCanyon) >= FT(1e-6)
        @warn "EBCanyon is not 0. Please check longwave_absorbed_no_tree"
    end

    if any(abs(getfield(LWREB_nT, field)) ≥ 1e-6 for field in propertynames(LWREB_nT))
        @warn "At least one LWREB_nT field is not 0. Please check longwave_absorbed_no_tree"
    end

    return LWRin_nT, LWRout_nT, LWRabs_nT, LWREB_nT
end

"""
    longwave_absorbed_with_trees(
            h_can::FT, w_can::FT, r_tree::FT, LWR::FT,
            fgveg::FT, fgbare::FT, fgimp::FT,
            ew::FT, et::FT, egveg::FT, egbare::FT, egimp::FT,
            Tgimp::FT, Tgbare::FT, Tgveg::FT,
            Twsun::FT, Twshade::FT, Ttree::FT,
            view_factor::ViewFactor{FT}
    ) where {FT<:AbstractFloat}

Calculate longwave radiation exchange in an urban canyon with trees.

# Arguments
- `h_can`: Building height [-]
- `w_can`: Ground width [-]
- `r_tree`: Size of the tree crown, crown radius [-]
- `LWR`: Atmospheric longwave radiation [W/m²]
- `fgveg`: Partitioning ground vegetation [-]
- `fgbare`: Partitioning ground bare [-]
- `fgimp`: Partitioning ground impervious [-]
- `ew`: Wall emissivity [-]
- `et`: Tree emissivity [-]
- `egveg`: Ground vegetation emissivity [-]
- `egbare`: Ground bare emissivity [-]
- `egimp`: Ground impervious emissivity [-]
- `Tgimp`: Temperature of ground impervious [K]
- `Tgbare`: Temperature of ground bare [K]
- `Tgveg`: Temperature of ground vegetated [K]
- `Twsun`: Temperature of wall sun [K]
- `Twshade`: Temperature of wall shade [K]
- `Ttree`: Temperature of tree [K]
- `view_factor`: View factors for radiation exchange

# Returns
- Tuple of (LWRin_T, LWRout_T, LWRabs_T, LWREB_T) containing longwave radiation components
"""
function longwave_absorbed_with_trees(
    h_can::FT,
    w_can::FT,
    r_tree::FT,
    LWR::FT,
    fgveg::FT,
    fgbare::FT,
    fgimp::FT,
    ew::FT,
    et::FT,
    egveg::FT,
    egbare::FT,
    egimp::FT,
    Tgimp::FT,
    Tgbare::FT,
    Tgveg::FT,
    Twsun::FT,
    Twshade::FT,
    Ttree::FT,
    view_factor::RayTracing.ViewFactor{FT},
) where {FT<:AbstractFloat}

    # Stefan-Boltzmann constant [W⋅m⁻²⋅K⁻⁴]
    bolzm = FT(5.67e-8)

    # Surface areas
    A_g = w_can
    A_w = h_can
    A_t = 4π * r_tree  # Area of two trees (2 * 2πr)

    # Check if view factors add up to 1
    SVF = [
        view_factor.F_gs_T + view_factor.F_gt_T + 2view_factor.F_gw_T,
        view_factor.F_ww_T + view_factor.F_wt_T + view_factor.F_wg_T + view_factor.F_ws_T,
        view_factor.F_sg_T + 2view_factor.F_sw_T + view_factor.F_st_T,
        view_factor.F_ts_T + 2view_factor.F_tw_T + view_factor.F_tt_T + view_factor.F_tg_T,
    ]

    all(x -> FT(0.999) ≤ x ≤ FT(1.001), SVF) || @warn "View factors do not add up to 1"

    # Surface indicators
    Cimp = fgimp > 0
    Cbare = fgbare > 0
    Cveg = fgveg > 0

    # View factor matrix for infinite reflections
    Tij = [
        1 0 0 -(1-egveg)*view_factor.F_gw_T*Cveg -(1-egveg)*view_factor.F_gw_T*Cveg -(1-egveg)*view_factor.F_gt_T*Cveg -(1-egveg)*view_factor.F_gs_T*Cveg;
        0 1 0 -(1-egbare)*view_factor.F_gw_T*Cbare -(1-egbare)*view_factor.F_gw_T*Cbare -(1-egbare)*view_factor.F_gt_T*Cbare -(1-egbare)*view_factor.F_gs_T*Cbare;
        0 0 1 -(1-egimp)*view_factor.F_gw_T*Cimp -(1-egimp)*view_factor.F_gw_T*Cimp -(1-egimp)*view_factor.F_gt_T*Cimp -(1-egimp)*view_factor.F_gs_T*Cimp;
        -(1-ew)*view_factor.F_wg_T*fgveg*Cveg -(1-ew)*view_factor.F_wg_T*fgbare*Cbare -(1-ew)*view_factor.F_wg_T*fgimp*Cimp 1 -(1-ew)*view_factor.F_ww_T -(1-ew)*view_factor.F_wt_T -(1-ew)*view_factor.F_ws_T;
        -(1-ew)*view_factor.F_wg_T*fgveg*Cveg -(1-ew)*view_factor.F_wg_T*fgbare*Cbare -(1-ew)*view_factor.F_wg_T*fgimp*Cimp -(1-ew)*view_factor.F_ww_T 1 -(1-ew)*view_factor.F_wt_T -(1-ew)*view_factor.F_ws_T;
        -(1-et)*view_factor.F_tg_T*fgveg*Cveg -(1-et)*view_factor.F_tg_T*fgbare*Cbare -(1-et)*view_factor.F_tg_T*fgimp*Cimp -(1-et)*view_factor.F_tw_T -(1-et)*view_factor.F_tw_T 1-(1-et)*view_factor.F_tt_T -(1-et)*view_factor.F_ts_T;
        0 0 0 0 0 0 1
    ]

    # Emitted radiation per surface
    Omega_i = [
        egveg*bolzm*(Tgveg^4)*Cveg;
        egbare*bolzm*(Tgbare^4)*Cbare;
        egimp*bolzm*(Tgimp^4)*Cimp;
        ew*bolzm*(Twsun^4);
        ew*bolzm*(Twshade^4);
        et*bolzm*(Ttree^4);
        LWR
    ]

    # Solve for outgoing radiation
    B_i = Tij \ Omega_i

    B_i[7] ≈ LWR || @warn "LWR mismatch after matrix inversion"

    # Matrix for incoming radiation
    Tij2 = [
        0 0 0 view_factor.F_gw_T*Cveg view_factor.F_gw_T*Cveg view_factor.F_gt_T*Cveg view_factor.F_gs_T*Cveg;
        0 0 0 view_factor.F_gw_T*Cbare view_factor.F_gw_T*Cbare view_factor.F_gt_T*Cbare view_factor.F_gs_T*Cbare;
        0 0 0 view_factor.F_gw_T*Cimp view_factor.F_gw_T*Cimp view_factor.F_gt_T*Cimp view_factor.F_gs_T*Cimp;
        view_factor.F_wg_T*fgveg*Cveg view_factor.F_wg_T*fgbare*Cbare view_factor.F_wg_T*fgimp*Cimp 0 view_factor.F_ww_T view_factor.F_wt_T view_factor.F_ws_T;
        view_factor.F_wg_T*fgveg*Cveg view_factor.F_wg_T*fgbare*Cbare view_factor.F_wg_T*fgimp*Cimp view_factor.F_ww_T 0 view_factor.F_wt_T view_factor.F_ws_T;
        view_factor.F_tg_T*fgveg*Cveg view_factor.F_tg_T*fgbare*Cbare view_factor.F_tg_T*fgimp*Cimp view_factor.F_tw_T view_factor.F_tw_T view_factor.F_tt_T view_factor.F_ts_T;
        0 0 0 0 0 0 0
    ]

    # Calculate incoming and net radiation
    A_i = Tij2 * B_i
    e_i = [egveg; egbare; egimp; ew; ew; et; FT(1)]

    # Calculate net absorbed radiation
    Qnet_i = similar(B_i)
    for i in eachindex(e_i)
        if e_i[i] == 1
            Qnet_i[i] = A_i[i] - Omega_i[i]
        else
            Qnet_i[i] = (e_i[i]*B_i[i] - Omega_i[i])/(1-e_i[i])
        end
    end
    Qnet_i[7] = 0  # Sky has fixed emission

    # Create output structures
    LWRin_T = LongwaveRadiation{FT}(;
        GroundImp=A_i[3]*Cimp,
        GroundBare=A_i[2]*Cbare,
        GroundVeg=A_i[1]*Cveg,
        Tree=A_i[6],
        WallSun=A_i[4],
        WallShade=A_i[5],
        TotalGround=fgveg*A_i[1] + fgbare*A_i[2] + fgimp*A_i[3],
        TotalCanyon=A_i[1]*fgveg*A_g/A_g +
                    A_i[2]*fgbare*A_g/A_g +
                    A_i[3]*fgimp*A_g/A_g +
                    A_i[4]*A_w/A_g +
                    A_i[5]*A_w/A_g +
                    A_i[6]*A_t/A_g,
    )

    LWRout_T = LongwaveRadiation{FT}(;
        GroundImp=B_i[3]*Cimp,
        GroundBare=B_i[2]*Cbare,
        GroundVeg=B_i[1]*Cveg,
        Tree=B_i[6],
        WallSun=B_i[4],
        WallShade=B_i[5],
        TotalGround=fgveg*B_i[1] + fgbare*B_i[2] + fgimp*B_i[3],
        TotalCanyon=B_i[1]*fgveg*A_g/A_g +
                    B_i[2]*fgbare*A_g/A_g +
                    B_i[3]*fgimp*A_g/A_g +
                    B_i[4]*A_w/A_g +
                    B_i[5]*A_w/A_g +
                    B_i[6]*A_t/A_g,
    )

    LWRabs_T = LongwaveRadiation{FT}(;
        GroundImp=Qnet_i[3]*Cimp,
        GroundBare=Qnet_i[2]*Cbare,
        GroundVeg=Qnet_i[1]*Cveg,
        Tree=Qnet_i[6],
        WallSun=Qnet_i[4],
        WallShade=Qnet_i[5],
        TotalGround=fgveg*Qnet_i[1] + fgbare*Qnet_i[2] + fgimp*Qnet_i[3],
        TotalCanyon=Qnet_i[1]*fgveg*A_g/A_g +
                    Qnet_i[2]*fgbare*A_g/A_g +
                    Qnet_i[3]*fgimp*A_g/A_g +
                    Qnet_i[4]*A_w/A_g +
                    Qnet_i[5]*A_w/A_g +
                    Qnet_i[6]*A_t/A_g,
    )

    LWREB_T = LWRin_T - LWRout_T - LWRabs_T

    # Energy balance check variables
    TotalLWRref_to_atm =
        B_i[1]*view_factor.F_sg_T*fgveg +
        B_i[2]*view_factor.F_sg_T*fgbare +
        B_i[3]*view_factor.F_sg_T*fgimp +
        B_i[4]*view_factor.F_sw_T +
        B_i[5]*view_factor.F_sw_T +
        B_i[6]*view_factor.F_st_T

    # Energy balance checks
    EBSurface = LWRin_T.TotalCanyon - LWRabs_T.TotalCanyon - LWRout_T.TotalCanyon
    EBCanyon = LWR - LWRabs_T.TotalCanyon - TotalLWRref_to_atm

    if abs(EBSurface) >= FT(1e-6)
        @warn "EBSurface is not 0. Please check longwave_absorbed_with_trees"
    end
    if abs(EBCanyon) >= FT(1e-6)
        error("EBCanyon is not 0. Please check longwave_absorbed_with_trees")
    end

    if any(abs(getfield(LWREB_T, field)) ≥ 1e-6 for field in propertynames(LWREB_T))
        @warn "At least one LWREB_T field is not 0. Please check longwave_absorbed_with_trees"
    end

    return LWRin_T, LWRout_T, LWRabs_T, LWREB_T
end
