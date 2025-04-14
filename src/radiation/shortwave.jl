"""
    direct_shortwave_trees(h_can::FT, d_tree::FT, h_tree::FT, r_tree::FT, theta_Z::FT, theta_n::FT, SWR_dir::FT) where {FT<:AbstractFloat}

Calculate direct shortwave radiation received by two trees in an urban canyon.

# Arguments
- `h_can`: normalized building height [-]
- `d_tree`: location of trees in the canyon, tree-wall distance [-]
- `h_tree`: height of trees, vertical level at the crown center [-]
- `r_tree`: size of the tree crown, crown radius [-]
- `theta_Z`: solar zenith angle [rad]
- `theta_n`: difference between solar azimuth angle and canyon orientation [rad]
- `SWR_dir`: direct shortwave radiation W/m^2 of horizontal surfaces [W/m^2]

# Returns
- `SWR_tree1`: direct shortwave radiation received by tree 1 per m^2 tree surface [W/m^2]
- `SWR_tree2`: direct shortwave radiation received by tree 2 per m^2 tree surface [W/m^2]
"""
function direct_shortwave_trees(
    h_can::FT, d_tree::FT, h_tree::FT, r_tree::FT, theta_Z::FT, theta_n::FT, SWR_dir::FT
) where {FT<:AbstractFloat}
    # Correction for infeasible tree height and radius length
    if 2*r_tree >= h_can
        r_tree = h_can/2 - FT(0.000001)
        @warn "tree diameter is bigger than canyon height and is set to the canyon height"
    end
    if h_tree + r_tree >= h_can
        h_tree = h_can - r_tree - FT(0.000001)
        @warn "tree height is bigger than canyon height and is set to the canyon height"
    end

    Xsi = tan(theta_Z) * abs(sin(theta_n))

    denominator = (h_can - h_tree)^2 - r_tree^2
    tan_theta1 =
        (
            (1-d_tree)*(h_can-h_tree) +
            r_tree*sqrt((1-d_tree)^2 + (h_can-h_tree)^2 - r_tree^2)
        )/denominator
    tan_theta2 =
        (
            (1-d_tree)*(h_can-h_tree) -
            r_tree*sqrt((1-d_tree)^2 + (h_can-h_tree)^2 - r_tree^2)
        )/denominator
    tan_theta3 =
        (
            d_tree*(h_can-h_tree) + r_tree*sqrt(d_tree^2 + (h_can-h_tree)^2 - r_tree^2)
        )/denominator
    tan_theta4 =
        (
            d_tree*(h_can-h_tree) - r_tree*sqrt(d_tree^2 + (h_can-h_tree)^2 - r_tree^2)
        )/denominator

    # Helper function to calculate SWR for a tree
    function calc_tree_swr(Xsi::FT, tan_upper::FT, tan_lower::FT, offset::FT)
        if Xsi >= tan_upper
            return FT(0)
        elseif Xsi < tan_upper && Xsi >= tan_lower
            return SWR_dir *
                   (r_tree*sqrt(1+Xsi^2) + offset - (h_can-h_tree)*Xsi)/(2*π*r_tree)
        elseif Xsi < tan_lower
            return SWR_dir * (2*r_tree*sqrt(1+Xsi^2))/(2*π*r_tree)
        else
            return FT(0)
        end
    end

    # Calculate SWR for both trees
    SWR_tree1 = calc_tree_swr(Xsi, tan_theta1, tan_theta2, 1-d_tree)
    SWR_tree2 = calc_tree_swr(Xsi, tan_theta3, tan_theta4, d_tree)

    return SWR_tree1, SWR_tree2
end

"""
    direct_shortwave_surfaces(
        h_can::FT, w_can::FT, d_tree::FT, h_tree::FT, r_tree::FT,
        theta_Z::FT, theta_n::FT, SWR_dir::FT, LAIt::FT, trees::Bool,
        ParVegTree::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate direct shortwave radiation received by different urban canyon surfaces.

# Arguments
- `h_can`: normalized building height [-]
- `w_can`: normalized street width [-]
- `d_tree`: location of trees in the canyon, tree-wall distance [-]
- `h_tree`: height of trees, vertical level at the crown center [-]
- `r_tree`: size of the tree crown, crown radius [-]
- `theta_Z`: solar zenith angle [rad]
- `theta_n`: difference between solar azimuth angle and canyon orientation [rad]
- `SWR_dir`: direct shortwave radiation W/m^2 of horizontal surfaces [W/m^2]
- `LAIt`: leaf area index of the tree [-]
- `trees`: boolean indicating if trees are present
- `ParVegTree`: named tuple containing vegetation parameters, including Kopt

# Returns
- `SWRdir_g`: direct shortwave radiation received by the ground [W/m^2]
- `SWRdir_wsun`: direct shortwave radiation received by the sunlit wall [W/m^2]
- `SWRdir_wshd`: direct shortwave radiation received by the shaded wall [W/m^2]
- `SWRdir_t`: direct shortwave radiation received by the tree [W/m^2]
"""
function direct_shortwave_surfaces(
    h_can::FT,
    w_can::FT,
    d_tree::FT,
    h_tree::FT,
    r_tree::FT,
    theta_Z::FT,
    theta_n::FT,
    SWR_dir::FT,
    LAIt::FT,
    trees::Bool,
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
) where {FT<:AbstractFloat}

    # Calculate tree-related parameters
    if !trees
        tau = zero(FT)
        SWRdir_t = zero(FT)
        X_Shadow, X_Tree, n_Shadow, n_Tree = shadow_length_no_tree(
            h_can, w_can, theta_Z, theta_n
        )
    else
        SWR_tree1, SWR_tree2 = direct_shortwave_trees(
            h_can, d_tree, h_tree, r_tree, theta_Z, theta_n, SWR_dir
        )
        tau = exp(-ParVegTree.Kopt * LAIt)
        SWRdir_t = (1 - tau) * (SWR_tree1 + SWR_tree2) / 2
        X_Shadow, X_Tree, n_Shadow, n_Tree = shadow_length_with_trees(
            h_can, w_can, d_tree, h_tree, r_tree, theta_Z, theta_n
        )
    end

    # Calculate radiation components
    Xsi = tan(theta_Z) * abs(sin(theta_n))
    SWRdir_g = SWR_dir * (1 - X_Shadow + tau * X_Tree)
    SWRdir_wsun = SWR_dir * Xsi * (1 - n_Shadow + tau * n_Tree)
    SWRdir_wshd = zero(FT)

    # Energy conservation check and correction
    A_g = w_can
    A_w = h_can
    A_t = 2 * π * r_tree

    check_total =
        A_g/A_g * SWRdir_g +
        A_w/A_g * SWRdir_wsun +
        A_w/A_g * SWRdir_wshd +
        2 * A_t/A_g * SWRdir_t

    if abs(check_total - SWR_dir) > 1e-10
        delta = check_total - SWR_dir
        SWRdir_t -= delta * A_g / (2 * A_t)
    end

    return SWRdir_g, SWRdir_wsun, SWRdir_wshd, SWRdir_t
end

"""
    shortwave_absorbed_no_trees(
        h_can::FT, w_can::FT,
        fgveg::FT, fgbare::FT, fgimp::FT,
        awraw::FT, agveg::FT, agbare::FT, agimp::FT,
        SWR_dir::FT, SWR_diff::FT,
        theta_z::FT, theta_n::FT,
        view_factor::ViewFactor{FT},
        ParVegTree::VegetationParams{FT},
        ParWindows::WindowParams{FT},
        bem_enabled::Bool
    ) where {FT<:AbstractFloat}

Calculate shortwave radiation exchange in an urban canyon without trees.
"""
function shortwave_absorbed_no_trees(
    h_can::FT,
    w_can::FT,
    fgveg::FT,
    fgbare::FT,
    fgimp::FT,
    awraw::FT,
    agveg::FT,
    agbare::FT,
    agimp::FT,
    SWR_dir::FT,
    SWR_diff::FT,
    theta_z::FT,
    theta_n::FT,
    view_factor::RayTracing.ViewFactor{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    bem_enabled::Bool,
) where {FT<:AbstractFloat}
    A_s = w_can
    A_g = h_can
    A_w = h_can

    # Check view factors sum to 1
    svf = [
        view_factor.F_gs_nT + 2view_factor.F_gw_nT,
        view_factor.F_ww_nT + view_factor.F_wg_nT + view_factor.F_ws_nT,
        view_factor.F_sg_nT + 2view_factor.F_sw_nT,
    ]

    all(x -> FT(0.999) ≤ x ≤ FT(1.001), svf) || @warn "View factors do not add up to 1"

    # Get direct shortwave radiation
    SWRdir_ground, SWRdir_wallsun, SWRdir_wallshade, _ = direct_shortwave_surfaces(
        h_can, w_can, NaN, NaN, NaN, theta_z, theta_n, SWR_dir, NaN, false, ParVegTree
    )

    # Surface indicators
    Cimp = fgimp > 0
    Cbare = fgbare > 0
    Cveg = fgveg > 0

    # Calculate wall albedo accounting for windows if BEM is enabled
    aw = if bem_enabled
        awraw*(1-ParWindows.GlazingRatio) + ParWindows.SolarAlbedo*ParWindows.GlazingRatio
    else
        awraw
    end

    # Albedos vector [veg, bare, imp, sun wall, shade wall, sky]
    ai = [agveg, agbare, agimp, aw, aw, FT(0)]

    # View factor matrices for infinite reflections
    Tij = [
        1 0 0 -agveg*view_factor.F_gw_nT*Cveg -agveg*view_factor.F_gw_nT*Cveg -agveg*view_factor.F_gs_nT*Cveg;
        0 1 0 -agbare*view_factor.F_gw_nT*Cbare -agbare*view_factor.F_gw_nT*Cbare -agbare*view_factor.F_gs_nT*Cbare;
        0 0 1 -agimp*view_factor.F_gw_nT*Cimp -agimp*view_factor.F_gw_nT*Cimp -agimp*view_factor.F_gs_nT*Cimp;
        -aw*view_factor.F_wg_nT*fgveg*Cveg -aw*view_factor.F_wg_nT*fgbare*Cbare -aw*view_factor.F_wg_nT*fgimp*Cimp 1 -aw*view_factor.F_ww_nT -aw*view_factor.F_ws_nT;
        -aw*view_factor.F_wg_nT*fgveg*Cveg -aw*view_factor.F_wg_nT*fgbare*Cbare -aw*view_factor.F_wg_nT*fgimp*Cimp -aw*view_factor.F_ww_nT 1 -aw*view_factor.F_ws_nT;
        0 0 0 0 0 1
    ]

    # Incoming shortwave radiation from sky
    Omega_i = [
        agveg*SWRdir_ground*Cveg;
        agbare*SWRdir_ground*Cbare;
        agimp*SWRdir_ground*Cimp;
        aw*SWRdir_wallsun;
        FT(0);
        SWR_diff
    ]

    # Solve for outgoing radiation
    B_i = Tij \ Omega_i

    # Check sky radiation preserved
    B_i[6] ≈ SWR_diff ||
        @warn "Incoming longwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal"

    # Matrix for incoming radiation
    Tij2 = [
        0 0 0 view_factor.F_gw_nT*Cveg view_factor.F_gw_nT*Cveg view_factor.F_gs_nT*Cveg;
        0 0 0 view_factor.F_gw_nT*Cbare view_factor.F_gw_nT*Cbare view_factor.F_gs_nT*Cbare;
        0 0 0 view_factor.F_gw_nT*Cimp view_factor.F_gw_nT*Cimp view_factor.F_gs_nT*Cimp;
        view_factor.F_wg_nT*fgveg*Cveg view_factor.F_wg_nT*fgbare*Cbare view_factor.F_wg_nT*fgimp*Cimp 0 view_factor.F_ww_nT view_factor.F_ws_nT;
        view_factor.F_wg_nT*fgveg*Cveg view_factor.F_wg_nT*fgbare*Cbare view_factor.F_wg_nT*fgimp*Cimp view_factor.F_ww_nT 0 view_factor.F_ws_nT;
        0 0 0 0 0 0
    ]

    SWRdir_i = [
        SWRdir_ground*Cveg;
        SWRdir_ground*Cbare;
        SWRdir_ground*Cimp;
        SWRdir_wallsun;
        FT(0);
        FT(0)
    ]

    # Calculate incoming radiation
    A_i1 = Tij2*B_i + SWRdir_i
    A_i = B_i ./ ai
    A_i[ai .== 0] .= A_i1[ai .== 0]
    A_i[6] = 0  # Sky has fixed emission

    # Calculate net absorbed radiation
    Q_net = A_i .- B_i

    # Energy balance checks
    SWRin_atm = SWR_dir + SWR_diff
    TotalSWRref_to_atm =
        B_i[1]*view_factor.F_sg_nT*fgveg +
        B_i[2]*view_factor.F_sg_nT*fgbare +
        B_i[3]*view_factor.F_sg_nT*fgimp +
        B_i[4]*view_factor.F_sw_nT +
        B_i[5]*view_factor.F_sw_nT

    canyon_albedo = TotalSWRref_to_atm/SWRin_atm

    # Create SWR output structure
    # Checked and correct
    SWRin_nT = LongwaveRadiation{FT}(;
        GroundImp=A_i[3]*Cimp,
        GroundBare=A_i[2]*Cbare,
        GroundVeg=A_i[1]*Cveg,
        Tree=FT(0),
        WallSun=A_i[4],
        WallShade=A_i[5],
        TotalGround=fgveg*A_i[1] + fgbare*A_i[2] + fgimp*A_i[3],
        TotalCanyon=A_i[1]*fgveg +
                    A_i[2]*fgbare +
                    A_i[3]*fgimp +
                    A_i[4]*h_can/w_can +
                    A_i[5]*h_can/w_can,
    )

    # Create additional radiation components using the same structure
    SWRout_nT = LongwaveRadiation{FT}(;
        GroundImp=B_i[3]*Cimp,
        GroundBare=B_i[2]*Cbare,
        GroundVeg=B_i[1]*Cveg,
        Tree=FT(0),
        WallSun=B_i[4],
        WallShade=B_i[5],
        TotalGround=fgveg*B_i[1] + fgbare*B_i[2] + fgimp*B_i[3],
        TotalCanyon=B_i[1]*fgveg +
                    B_i[2]*fgbare +
                    B_i[3]*fgimp +
                    B_i[4]*h_can/w_can +
                    B_i[5]*h_can/w_can,
    )

    SWRabs_nT = LongwaveRadiation{FT}(;
        GroundImp=Q_net[3]*Cimp,
        GroundBare=Q_net[2]*Cbare,
        GroundVeg=Q_net[1]*Cveg,
        Tree=FT(0),
        WallSun=Q_net[4],
        WallShade=Q_net[5],
        TotalGround=fgveg*Q_net[1] + fgbare*Q_net[2] + fgimp*Q_net[3],
        TotalCanyon=Q_net[1]*fgveg +
                    Q_net[2]*fgbare +
                    Q_net[3]*fgimp +
                    Q_net[4]*h_can/w_can +
                    Q_net[5]*h_can/w_can,
    )

    # Direct absorbed shortwave radiation
    SWRabsDir_nT = LongwaveRadiation{FT}(;
        GroundImp=(1-agimp)*SWRdir_ground*Cimp,
        GroundBare=(1-agbare)*SWRdir_ground*Cbare,
        GroundVeg=(1-agveg)*SWRdir_ground*Cveg,
        Tree=FT(0),
        WallSun=(1-aw)*SWRdir_wallsun,
        WallShade=(1-aw)*SWRdir_wallshade,
        TotalGround=fgveg*(1-agveg)*SWRdir_ground +
                    fgbare*(1-agbare)*SWRdir_ground +
                    fgimp*(1-agimp)*SWRdir_ground,
        TotalCanyon=fgveg*(1-agveg)*SWRdir_ground*A_g/A_g +
                    fgbare*(1-agbare)*SWRdir_ground*A_g/A_g +
                    fgimp*(1-agimp)*SWRdir_ground*A_g/A_g +
                    (1-aw)*SWRdir_wallsun*A_w/A_g +
                    (1-aw)*SWRdir_wallshade*A_w/A_g,
    )

    # Diffuse absorbed shortwave radiation
    SWRabsDiff_nT = LongwaveRadiation{FT}(;
        GroundImp=(SWRabs_nT.GroundImp - SWRabsDir_nT.GroundImp)*Cimp,
        GroundBare=(SWRabs_nT.GroundBare - SWRabsDir_nT.GroundBare)*Cbare,
        GroundVeg=(SWRabs_nT.GroundVeg - SWRabsDir_nT.GroundVeg)*Cveg,
        Tree=FT(0),
        WallSun=(SWRabs_nT.WallSun - SWRabsDir_nT.WallSun),
        WallShade=(SWRabs_nT.WallShade - SWRabsDir_nT.WallShade),
        TotalGround=(SWRabs_nT.TotalGround - SWRabsDir_nT.TotalGround),
        TotalCanyon=(SWRabs_nT.TotalCanyon - SWRabsDir_nT.TotalCanyon),
    )

    # Energy Balance of shortwave radiation
    SWREB_nT = SWRin_nT - SWRout_nT - SWRabs_nT

    # Energy balance checks
    for field in propertynames(SWREB_nT)
        val = getfield(SWREB_nT, field)
        if abs(val) >= FT(1e-6)
            @warn "$(String(field)) is not 0. Please check shortwave_absorbed_no_trees"
        end
    end

    return SWRin_nT,
    SWRout_nT, SWRabs_nT, SWRabsDir_nT, SWRabsDiff_nT, SWREB_nT,
    canyon_albedo
end

"""
    shortwave_absorbed_with_trees(
            h_can::FT, w_can::FT, h_tree::FT, r_tree::FT, d_tree::FT,
            fgveg::FT, fgbare::FT, fgimp::FT,
            awraw::FT, agveg::FT, agbare::FT, agimp::FT, at::FT,
            LAIt::FT, SWR_dir::FT, SWR_diff::FT,
            theta_Z::FT, theta_n::FT,
            view_factor::ViewFactor{FT},
            ParVegTree::VegetationParams{FT},
            ParWindows::WindowParams{FT},
            bem_enabled::Bool
    ) where {FT<:AbstractFloat}

Calculate shortwave radiation exchange in an urban canyon with trees.

# Returns
Tuple of (SWRin_T, SWRout_T, SWRabs_T, SWRabsDir_T, SWRabsDiff_T, SWREB_T, albedo_canyon)
"""
function shortwave_absorbed_with_trees(
    h_can::FT,
    w_can::FT,
    h_tree::FT,
    r_tree::FT,
    d_tree::FT,
    fgveg::FT,
    fgbare::FT,
    fgimp::FT,
    awraw::FT,
    agveg::FT,
    agbare::FT,
    agimp::FT,
    at::FT,
    LAIt::FT,
    SWR_dir::FT,
    SWR_diff::FT,
    theta_Z::FT,
    theta_n::FT,
    view_factor::RayTracing.ViewFactor{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    bem_enabled::Bool,
) where {FT<:AbstractFloat}

    # Normalized surface areas
    A_s = w_can
    A_g = w_can
    A_w = h_can
    A_t = 4π * r_tree  # Area of two trees (2 * 2πr)

    # Calculate direct shortwave radiation components
    SWRdir_ground, SWRdir_wallsun, SWRdir_wallshade, SWRdir_tree = direct_shortwave_surfaces(
        h_can,
        w_can,
        d_tree,
        h_tree,
        r_tree,
        theta_Z,
        theta_n,
        SWR_dir,
        LAIt,
        true,
        ParVegTree,
    )

    # Surface indicators
    Cimp = fgimp > 0
    Cbare = fgbare > 0
    Cveg = fgveg > 0

    # Calculate wall albedo accounting for windows if BEM is enabled
    aw = if bem_enabled
        awraw * (1-ParWindows.GlazingRatio) + ParWindows.SolarAlbedo * ParWindows.GlazingRatio
    else
        awraw
    end

    # Albedos vector [veg, bare, imp, sun wall, shade wall, trees, sky]
    ai = [agveg, agbare, agimp, aw, aw, at, FT(0)]

    # View factor matrix for infinite reflections
    Tij = [
        1 0 0 -agveg*view_factor.F_gw_T*Cveg -agveg*view_factor.F_gw_T*Cveg -agveg*view_factor.F_gt_T*Cveg -agveg*view_factor.F_gs_T*Cveg;
        0 1 0 -agbare*view_factor.F_gw_T*Cbare -agbare*view_factor.F_gw_T*Cbare -agbare*view_factor.F_gt_T*Cbare -agbare*view_factor.F_gs_T*Cbare;
        0 0 1 -agimp*view_factor.F_gw_T*Cimp -agimp*view_factor.F_gw_T*Cimp -agimp*view_factor.F_gt_T*Cimp -agimp*view_factor.F_gs_T*Cimp;
        -aw*view_factor.F_wg_T*fgveg*Cveg -aw*view_factor.F_wg_T*fgbare*Cbare -aw*view_factor.F_wg_T*fgimp*Cimp 1 -aw*view_factor.F_ww_T -aw*view_factor.F_wt_T -aw*view_factor.F_ws_T;
        -aw*view_factor.F_wg_T*fgveg*Cveg -aw*view_factor.F_wg_T*fgbare*Cbare -aw*view_factor.F_wg_T*fgimp*Cimp -aw*view_factor.F_ww_T 1 -aw*view_factor.F_wt_T -aw*view_factor.F_ws_T;
        -at*view_factor.F_tg_T*fgveg*Cveg -at*view_factor.F_tg_T*fgbare*Cbare -at*view_factor.F_tg_T*fgimp*Cimp -at*view_factor.F_tw_T -at*view_factor.F_tw_T 1-at*view_factor.F_tt_T -at*view_factor.F_ts_T;
        0 0 0 0 0 0 1
    ]

    # Incoming shortwave radiation from sky
    Omega_i = [
        agveg*SWRdir_ground*Cveg;
        agbare*SWRdir_ground*Cbare;
        agimp*SWRdir_ground*Cimp;
        aw*SWRdir_wallsun;
        aw*zero(FT);
        at*SWRdir_tree;
        SWR_diff
    ]

    # Solve for outgoing radiation
    B_i = Tij \ Omega_i

    B_i[7] ≈ SWR_diff || @warn "SWR_diff mismatch after matrix inversion"

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

    SWRdir_i = [
        SWRdir_ground*Cveg;
        SWRdir_ground*Cbare;
        SWRdir_ground*Cimp;
        SWRdir_wallsun;
        zero(FT);
        SWRdir_tree;
        zero(FT)
    ]

    # Calculate incoming radiation
    A_i1 = Tij2*B_i + SWRdir_i
    A_i = B_i ./ ai
    A_i[ai .== 0] .= A_i1[ai .== 0]
    A_i[7] = 0  # Sky has fixed emission

    # Calculate net absorbed radiation
    Q_net = A_i .- B_i

    # Create output structures
    SWRin_T = create_longwave_radiation(
        A_i, fgveg, fgbare, fgimp, h_can, w_can, A_t, Cimp, Cbare, Cveg
    )
    SWRout_T = create_longwave_radiation(
        B_i, fgveg, fgbare, fgimp, h_can, w_can, A_t, Cimp, Cbare, Cveg
    )
    SWRabs_T = create_longwave_radiation(
        Q_net, fgveg, fgbare, fgimp, h_can, w_can, A_t, Cimp, Cbare, Cveg
    )

    # Direct absorbed shortwave radiation components
    SWRabsDir_T = create_longwave_radiation(
        [
            (1-agveg)*SWRdir_ground*Cveg;
            (1-agbare)*SWRdir_ground*Cbare;
            (1-agimp)*SWRdir_ground*Cimp;
            (1-aw)*SWRdir_wallsun;
            (1-aw)*SWRdir_wallshade;
            (1-at)*SWRdir_tree;
            zero(FT)
        ],
        fgveg,
        fgbare,
        fgimp,
        h_can,
        w_can,
        A_t,
        Cimp,
        Cbare,
        Cveg,
    )

    # Diffuse absorbed shortwave radiation
    SWRabsDiff_T = SWRabs_T - SWRabsDir_T

    # Energy balance
    SWREB_T = SWRin_T - SWRout_T - SWRabs_T

    # Calculate canyon albedo
    TotalSWRref_to_atm =
        B_i[1]*view_factor.F_sg_T*fgveg +
        B_i[2]*view_factor.F_sg_T*fgbare +
        B_i[3]*view_factor.F_sg_T*fgimp +
        B_i[4]*view_factor.F_sw_T +
        B_i[5]*view_factor.F_sw_T +
        B_i[6]*view_factor.F_st_T

    albedo_canyon = TotalSWRref_to_atm/(SWR_dir + SWR_diff)

    # Energy balance checks
    if any(abs(getfield(SWREB_T, field)) ≥ FT(1e-6) for field in propertynames(SWREB_T))
        error("At least one SWREB_T field is not 0")
    end

    return SWRin_T, SWRout_T, SWRabs_T, SWRabsDir_T, SWRabsDiff_T, SWREB_T, albedo_canyon
end
