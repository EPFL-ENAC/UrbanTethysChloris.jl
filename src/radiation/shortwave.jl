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
