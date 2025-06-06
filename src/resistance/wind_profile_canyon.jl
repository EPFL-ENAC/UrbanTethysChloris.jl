"""
    wind_profile_canyon(
        Hcan::FT,
        Htree::FT,
        R_tree::FT,
        Wcan::FT,
        Wroof::FT,
        Kopt::FT,
        LAI_t::FT,
        Zatm::FT,
        uatm::FT,
        Zp::FT,
        trees::Bool,
        Zref_und::FT,
        zom_und::FT,
        Hcan_max::FT,
        Hcan_std::FT
    ) where {FT<:AbstractFloat}

Calculate urban canyon wind profile parameters.

# Arguments
- `Hcan`: Canyon height [m]
- `Htree`: Tree height [m]
- `R_tree`: Tree radius [m]
- `Wcan`: Canyon width [m]
- `Wroof`: Roof width [m]
- `Kopt`: Optical transmission factor [-]
- `LAI_t`: Tree leaf area index [-]
- `Zatm`: Atmospheric reference height [m]
- `uatm`: Wind speed at atmospheric reference height [m/s]
- `Zp`: Height of interest within canyon [m]
- `trees`: Presence of trees
- `Zref_und`: Reference height within canyon [m]
- `zom_und`: Roughness length of ground surface [m]
- `Hcan_max`: Maximum canyon height [m]
- `Hcan_std`: Standard deviation of canyon height [m]

# Returns
- `dcan`: Urban displacement height including trees [m]
- `zomcan`: Urban momentum roughness height including trees [m]
- `u_Hcan_max`: Wind speed at canyon height [m/s]
- `u_Zp`: Wind speed within canyon at height Zp [m/s]
- `w_Zp`: Vertical wind speed within canyon [m/s]
- `alpha`: Canyon attenuation coefficient [-]
- `RoughnessParameter`: Type of roughness parameterization used
"""
function wind_profile_canyon(
    Hcan::FT,
    Htree::FT,
    R_tree::FT,
    Wcan::FT,
    Wroof::FT,
    Kopt::FT,
    LAI_t::FT,
    Zatm::FT,
    uatm::FT,
    Zp::FT,
    trees::Bool,
    Zref_und::FT,
    zom_und::FT,
    Hcan_max::FT,
    Hcan_std::FT,
) where {FT<:AbstractFloat}
    # Adjust for no trees case
    if R_tree == 0
        trees = false
    end

    if !trees
        Htree = zero(FT)
        R_tree = zero(FT)
        LAI_t = zero(FT)
    end

    # Parameters McDonald 1998
    a = FT(4.43)    # Best fit for staggered arrays
    k = FT(0.4)     # Von Karman constant
    b = FT(1.0)     # No drag correction factors
    CDb = FT(1.2)   # Nominal drag for cubical obstacles

    # Plan and frontal area fractions
    Ap_build = Wroof
    Ap_tree = 4 * R_tree
    Ap_urb = Wcan + Wroof

    TreeCrownHeight = min(FT(3), FT(0.5) * Hcan)

    Af_build_s = Hcan * Wroof
    Af_veg_s = 4 * R_tree * TreeCrownHeight
    Ap_urbArea = (Wcan + Wroof) * Wroof

    # Tree canopy transmittance
    P2D = exp(-Kopt * LAI_t)
    P3D = P2D^FT(0.40)
    Pv = (-FT(1.251) * P3D^2 + FT(0.489) * P3D + FT(0.803)) / CDb

    # Structural parameters
    Lp_tot = (Ap_build + (1 - P3D) * Ap_tree) / Ap_urb
    H_tot =
        (Hcan * Ap_build + (Htree + R_tree) * (1 - P3D) * Ap_tree) /
        (Ap_build + (1 - P3D) * Ap_tree)

    # MacDonald calculations
    dcan_MacD = H_tot * (1 + a^(-Lp_tot) * (Lp_tot - 1))

    Af_build = Af_build_s
    Af_veg = Af_veg_s

    zomcan_MacD =
        H_tot *
        (1 - dcan_MacD/H_tot) *
        exp(
            -(
                1/k^2 *
                FT(0.5) *
                b *
                CDb *
                (1 - dcan_MacD/H_tot) *
                (Af_build + Pv * Af_veg) / Ap_urbArea
            )^(-FT(0.5)),
        )

    # Kanda corrections if height variability data available
    if !isnan(Hcan_max) && !isnan(Hcan_std)
        a0 = FT(1.29)
        b0 = FT(0.36)
        c0 = FT(-0.17)

        X = (Hcan_std + H_tot) / Hcan_max
        dcan_Kand = (c0 * X^2 + (a0 * Lp_tot^b0 - c0) * X) * Hcan_max

        a1 = FT(0.71)
        b1 = FT(20.21)
        c1 = FT(-0.77)

        Y = Lp_tot * Hcan_std / H_tot
        zomcan_Kand = (b1 * Y^2 + c1 * Y + a1) * zomcan_MacD

        if (dcan_Kand + zomcan_Kand) >= (Hcan_max - FT(0.5))
            dcan = dcan_MacD
            zomcan = zomcan_MacD
            Hcan_max = Hcan
            RoughnessParameter = :MacD
        else
            dcan = dcan_Kand
            zomcan = zomcan_Kand
            RoughnessParameter = :Kand
        end
    else
        dcan = dcan_MacD
        zomcan = zomcan_MacD
        Hcan_max = Hcan
        RoughnessParameter = :MacD
    end

    # Wind profile calculations
    us_atm = k * uatm / log((Zatm - dcan)/zomcan)
    u_Hcan_max = (us_atm/k) * log((Hcan_max - dcan)/zomcan)
    alpha = log(uatm/u_Hcan_max)/(Zatm/Hcan_max - 1)

    # Calculate wind speeds at different heights
    if Zp >= Hcan_max
        u_Zp = (us_atm/k) * log((Zp - dcan)/zomcan)
        w_Zp = zero(FT)
    elseif Zp <= Hcan_max && Zp >= Zref_und
        u_Zp = u_Hcan_max * exp(-alpha * (1 - Zp/Hcan_max))
        w_Zp = zero(FT)
    elseif Zp <= Zref_und && Zp >= zom_und
        uref_und = u_Hcan_max * exp(-alpha * (1 - Zref_und/Hcan_max))
        usref_und = k * uref_und / log(Zref_und/zom_und)
        u_Zp = (usref_und/k) * log(Zp/zom_und)
        w_Zp = zero(FT)
    else
        u_Zp = zero(FT)
        w_Zp = zero(FT)
        @info "wind speed calculation height higher than reference height or lower than roughness length"
    end

    return dcan, zomcan, u_Hcan_max, u_Zp, w_Zp, alpha, RoughnessParameter
end
