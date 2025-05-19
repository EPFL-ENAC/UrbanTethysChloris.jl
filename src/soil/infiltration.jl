"""
    infiltration(
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
        Pond::FT
    ) where {FT<:AbstractFloat}

Calculate infiltration rates based on soil parameters and conditions.

# Arguments
- `Osat`: Saturated water content [m³/m³]
- `Ohy`: Hygroscopic water content [m³/m³]
- `L`: Pore size distribution index [-]
- `alpVG`: van Genuchten α parameter [1/mm]
- `nVG`: van Genuchten n parameter [-]
- `Pe`: Air entry pressure [kPa]
- `Ks_Zs`: Saturated hydraulic conductivity [mm/h]
- `O33`: Water content at -33 kPa [m³/m³]
- `SPAR`: Soil parameterization choice:
    1. van Genuchten (1980)
    2. Saxton-Rawls (1986)
- `O`: Current soil water content [m³/m³]
- `Dz`: Distance from surface to half-layer [mm]
- `WIS`: Water incoming to soil layer [mm/h]
- `cosalp`: Cosine of slope angle [-]
- `Pond`: Ponding depth [mm]

# Returns
- `f`: Actual infiltration rate [mm/h]
- `fpot`: Potential infiltration rate [mm/h]
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
