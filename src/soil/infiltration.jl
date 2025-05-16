"""
    infiltration(Osat, Ohy, L, alpVG, nVG, lVG, Pe, Ks_Zs, O33, Ks_mac, Omac, alpVGM, nVGM, lVGM, Phy, s_SVG, bVG, SPAR, O, Dz, WIS, cosalp, Pond)

Calculate infiltration rates based on soil parameters and conditions.

# Arguments
- All input parameters should be of the same floating-point type
- `Osat`: Saturation moisture at 0 kPa
- `O`: Water Content Soil Active Layer for Infiltration
- `L`: Slope of logarithmic tension-moisture curve
- `Pe`: Tension at air entry (bubbling pressure) [kPa]
- `Ks_Zs`: Saturation conductivity [mm/h]
- `WIS`: Water Incoming to Soil Layer [mm/h]
- `Dz`: Distance from surface to half-layer [mm]
- `Pond`: Ponding depth [mm]

# Returns
- `f`: Infiltration rate [mm/h]
- `fpot`: Potential Infiltration rate [mm/h]
"""
function infiltration(
    Osat::FT,
    Ohy::FT,
    L::FT,
    alpVG::FT,
    nVG::FT,
    Pe::FT,
    Ks_Zs::FT,
    O33::FT,
    SPAR::Int,
    O::FT,
    Dz::FT,
    WIS::FT,
    cosalp::FT,
    Pond::FT,
) where {FT<:AbstractFloat}

    # Calculate PAE (suction negative) based on SPAR
    PAE = if SPAR == 2
        gw = FT(9810)
        -FT(1000) * FT(1000) * Pe / gw  # [mm]
    else
        zero(FT)  # Cases 1
    end

    # Handle ponding
    P0 = max(zero(FT), Pond)

    # Get conductivity and suction
    K, P = conductivity_suction(SPAR, Ks_Zs, Osat, Ohy, L, Pe, O33, alpVG, nVG, O)

    P = -P  # [mm] Water Potential
    Khalf = FT(0.5) * (K + Ks_Zs)

    # Calculate potential infiltration rate [mm/h]
    fpot = Khalf * (one(FT) * cosalp - (-PAE + (P - P0)) / Dz)

    # Calculate actual infiltration rate [mm/h]
    f = max(zero(FT), min(fpot, WIS))

    return f, fpot
end
