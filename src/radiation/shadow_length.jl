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
