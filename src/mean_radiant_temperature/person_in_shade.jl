"""
    person_in_shade(
        trees::Bool,
        h_can::FT,
        w_can::FT,
        d_tree::FT,
        h_tree::FT,
        r_tree::FT,
        theta_Z::FT,
        theta_n::FT,
        h_P::FT,
        x_P::FT,
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        Wcan::FT,
        TimeOfMaxSolAlt::FT,
        TimeHr::FT,
    ) where {FT<:AbstractFloat}

Determine if a person is in shade.

# Arguments
- `trees`: Trees present in canyon
- `h_can`: Normalized building height [-]
- `w_can`: Normalized street width [-]
- `d_tree`: Location of trees in canyon (tree-wall distance) [-]
- `h_tree`: Height of trees (vertical level at crown center) [-]
- `r_tree`: Size of tree crown (crown radius) [-]
- `theta_Z`: Solar zenith angle [rad]
- `theta_n`: Difference between solar azimuth angle and canyon orientation [rad]
- `h_P`: Vertical position of person H_p[m]/W_can[m] [-]
- `x_P`: Relative position within canyon (0=left edge, 1=right edge) [-]
- `ParVegTree`: Tree vegetation parameters
- `Wcan`: Canyon width [m]
- `TimeOfMaxSolAlt`: Time of maximum solar altitude [h]
- `TimeHr`: Current time [h]

# Returns
- `BooleanInSun::FT`: Factor indicating if person is in sun (0=shade, 1=sun, 0-1=partial)
"""
function person_in_shade(
    trees::Bool,
    h_can::FT,
    w_can::FT,
    d_tree::FT,
    h_tree::FT,
    r_tree::FT,
    theta_Z::FT,
    theta_n::FT,
    h_P::FT,
    x_P::FT,
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    Wcan::FT,
    TimeOfMaxSolAlt::FT,
    TimeHr::FT,
) where {FT<:AbstractFloat}

    # Normalize person height and position
    h_P = h_P / Wcan
    x_P = TimeHr <= TimeOfMaxSolAlt ? x_P / Wcan : (Wcan - x_P) / Wcan

    # Tree parameters
    Kopt_T = ParVegTree.Kopt
    LAI_T = ParVegTree.LAI

    # Correct for infeasible tree height and radius length
    if h_tree + r_tree >= h_can
        h_tree = h_can - r_tree - FT(0.000001)
        @warn "Tree height is bigger than canyon height and is set to the canyon height. The tree height is adjusted."
    end

    # Shift horizontal plane to point height
    h_can = h_can - h_P
    h_tree = h_tree - h_P

    # Calculate shadow locations
    Xsi = tan(theta_Z) * abs(sin(theta_n))
    secXsi = sqrt(1 + Xsi^2)

    # Shadow location by wall (origin (0,0) is lower left corner of canyon)
    x0 = max(zero(FT), w_can - h_can * Xsi)

    # Shadow location by Tree 1
    x1 = max(zero(FT), d_tree - h_tree * Xsi - r_tree * secXsi)
    x2 = max(zero(FT), d_tree - h_tree * Xsi + r_tree * secXsi)

    # Shadow location by Tree 2
    x3 = max(zero(FT), w_can - d_tree - h_tree * Xsi - r_tree * secXsi)
    x4 = max(zero(FT), w_can - d_tree - h_tree * Xsi + r_tree * secXsi)

    # Calculate if in shade or not
    if trees
        tau = exp(-Kopt_T * LAI_T)  # Shortwave radiation passing through trees
        if x_P > x0  # In full shade of wall
            BooleanInSun = zero(FT)
        elseif (x1 < x_P < x2)  # In shade of tree 1
            BooleanInSun = tau
        elseif (x3 < x_P < x4)  # In shade of tree 2
            BooleanInSun = tau
        else  # In sun
            BooleanInSun = one(FT)
        end
    else
        if x_P > x0  # In full shade of wall
            BooleanInSun = zero(FT)
        else  # In sun
            BooleanInSun = one(FT)
        end
    end

    return BooleanInSun
end
