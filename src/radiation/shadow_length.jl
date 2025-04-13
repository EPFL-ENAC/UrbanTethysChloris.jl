"""
    shadow_length_no_tree(h_can::FT, w_can::FT, theta_Z::FT, theta_n::FT) where {FT<:AbstractFloat}

Calculate shadow lengths in an urban canyon without trees.

# Arguments
- `h_can`: normalized building height [-]
- `w_can`: normalized street width [-]
- `theta_Z`: solar zenith angle [rad]
- `theta_n`: difference between solar azimuth angle and canyon orientation [rad]

# Returns
- `X_Shadow`: fraction of ground that is shaded [0-1]
- `X_tree`: fraction of ground that is shaded by tree [0-1]
- `n_Shadow`: fraction of wall that is shaded [0-1]
- `n_tree`: fraction of wall that is shaded by tree [0-1]
"""
function shadow_length_no_tree(
    h_can::FT, w_can::FT, theta_Z::FT, theta_n::FT
) where {FT<:AbstractFloat}
    # Compute shadow geometry
    Xsi = tan(theta_Z) * abs(sin(theta_n))

    # Shadow by the Wall
    X_Shadow = h_can * Xsi            # shadow cast on the ground by the wall
    n_Shadow = h_can - w_can/Xsi      # shadow cast on the opposite wall by the wall

    if abs(X_Shadow) < w_can
        n_Shadow = 0.0
    else
        X_Shadow = w_can
    end

    n_Shadow = min(h_can, n_Shadow)

    X_tree = 0.0
    n_tree = 0.0
    n_Shadow = n_Shadow/h_can

    return X_Shadow, X_tree, n_Shadow, n_tree
end

"""
    shadow_length_with_trees(h_can::FT, w_can::FT, d_tree::FT, h_tree::FT, r_tree::FT, theta_Z::FT, theta_n::FT) where {FT<:AbstractFloat}

Calculate shadow lengths in an urban canyon with trees.

# Arguments
- `h_can`: normalized building height [-]
- `w_can`: normalized street width [-]
- `d_tree`: location of trees in the canyon, tree-wall distance [-]
- `h_tree`: height of trees, vertical level at the crown center [-]
- `r_tree`: size of the tree crown, crown radius [-]
- `theta_Z`: solar zenith angle [rad]
- `theta_n`: difference between solar azimuth angle and canyon orientation [rad]

# Returns
- `X_Shadow`: fraction of ground that is shaded [0-1]
- `X_tree`: fraction of ground that is shaded by tree [0-1]
- `n_Shadow`: fraction of wall that is shaded [0-1]
- `n_tree`: fraction of wall that is shaded by tree [0-1]
"""
function shadow_length_with_trees(
    h_can::FT, w_can::FT, d_tree::FT, h_tree::FT, r_tree::FT, theta_Z::FT, theta_n::FT
) where {FT<:AbstractFloat}
    # Correction for infeasible tree height and radius length
    if 2*r_tree >= h_can
        r_tree = h_can/2 - 0.000001
        @warn "tree diameter is bigger than canyon height and is set to the canyon height. The radius is adjusted."
    end
    if h_tree + r_tree >= h_can
        h_tree = h_can - r_tree - 0.000001
        @warn "tree height is bigger than canyon height and is set to the canyon height. The tree height is adjusted."
    end

    # CALCULATION
    Xsi = tan(theta_Z) * abs(sin(theta_n))

    # Shadow by the Wall
    X_wall = h_can * Xsi
    n_wall = h_can - w_can/Xsi

    if abs(X_wall) < w_can
        n_wall = 0.0
    else
        X_wall = w_can
    end

    x0 = max(0.0, w_can - h_can*Xsi)
    y0 = max(0.0, h_can - w_can/Xsi)
    secXsi = sqrt(1 + Xsi^2)
    cosecXsi = sqrt(1 + 1/Xsi^2)

    # Shadow by the Tree 1
    x1 = max(0.0, d_tree - h_tree*Xsi - r_tree*secXsi)
    y1 = max(0.0, h_tree - (w_can-d_tree)/Xsi - r_tree*cosecXsi)
    x2 = max(0.0, d_tree - h_tree*Xsi + r_tree*secXsi)
    y2 = max(0.0, h_tree - (w_can-d_tree)/Xsi + r_tree*cosecXsi)

    X_Tree1 = x2 - x1

    # Shadow by the Tree 2
    x3 = max(0.0, w_can-d_tree - h_tree*Xsi - r_tree*secXsi)
    y3 = max(0.0, h_tree - d_tree/Xsi - r_tree*cosecXsi)
    x4 = max(0.0, w_can-d_tree - h_tree*Xsi + r_tree*secXsi)
    y4 = max(0.0, h_tree - d_tree/Xsi + r_tree*cosecXsi)

    X_Tree2 = x4 - x3
    n_Tree1 = y4 - y3
    n_Tree2 = y2 - y1

    # Total shadow length on the ground
    delta = max(0.0, x2-x0)

    if x0 < x4
        X_Shadow = w_can - min(x0, x3) + X_Tree1 - delta
        X_Tree = if x0 < x3
            X_Tree1 - delta
        else
            X_Tree1 + x0 - x3
        end
    else
        X_Shadow = X_wall + X_Tree1 + X_Tree2
        X_Tree = X_Tree1 + X_Tree2
    end

    # Total shadow length on the wall
    lowest_shaded = max(y0, y2)

    if y3 > lowest_shaded
        n_Shadow = n_Tree1 + lowest_shaded
        n_Tree = if y2 > n_wall
            n_Tree1 + y2 - n_wall
        elseif y2 <= n_wall
            n_Tree1
        elseif y1 > n_wall
            n_Tree1 + n_Tree2
        else
            n_Tree1
        end
    else
        n_Shadow = max(y0, y1, y2, y3, y4)
        n_Tree = if y4 > n_wall
            y4 - n_wall
        else
            0.0
        end
    end

    n_Tree = n_Tree/h_can
    n_Shadow = n_Shadow/h_can

    return X_Shadow, X_Tree, n_Shadow, n_Tree
end
