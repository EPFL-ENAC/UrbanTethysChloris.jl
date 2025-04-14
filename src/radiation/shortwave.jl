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
