"""
    soil_moisture_conductivity_update(
        V::Vector{FT},
        Pcla::FT,
        Psan::FT,
        Porg::FT,
        Kfc::FT,
        Phy::FT,
        SPAR::Int,
        Kbot::FT,
        CASE_ROOT_H::Int,
        CASE_ROOT_L::Int,
        ZR95_H::Vector{FT},
        ZR95_L::Vector{FT},
        ZR50_H::Vector{FT},
        ZR50_L::Vector{FT},
        ZRmax_H::Vector{FT},
        ZRmax_L::Vector{FT},
        Zs::Vector{FT},
        Rrootl_H::Vector{FT},
        Rrootl_L::Vector{FT},
        PsiL50_H::Vector{FT},
        PsiL50_L::Vector{FT},
        PsiX50_H::Vector{FT},
        PsiX50_L::Vector{FT}
    ) where {FT<:AbstractFloat}

Updates soil moisture content and calculates hydraulic properties based on current conditions.

# Arguments
- `V`: Soil water volume per unit area [mm]
- `Pcla`: Clay fraction in soil [-]
- `Psan`: Sand fraction in soil [-]
- `Porg`: Organic matter fraction in soil [-]
- `Kfc`: Hydraulic conductivity at field capacity [mm/h]
- `Phy`: Soil water potential at hygroscopic point [kPa]
- `SPAR`: Soil parameterization choice:
    1. Van-Genuchten (1980)
    2. Saxton-Rawls (1986)
- `Kbot`: Hydraulic conductivity at bottom boundary [mm/h]
- `CASE_ROOT_H`: Root distribution type for high vegetation [-]
- `CASE_ROOT_L`: Root distribution type for low vegetation [-]
- `ZR95_H`: 95th percentile root depth for high vegetation [mm]
- `ZR95_L`: 95th percentile root depth for low vegetation [mm]
- `ZR50_H`: 50th percentile root depth for high vegetation [mm]
- `ZR50_L`: 50th percentile root depth for low vegetation [mm]
- `ZRmax_H`: Maximum root depth for high vegetation [mm]
- `ZRmax_L`: Maximum root depth for low vegetation [mm]
- `Zs`: Soil layer depths [mm]
- `Rrootl_H`: Root length density for high vegetation [mm/mm³]
- `Rrootl_L`: Root length density for low vegetation [mm/mm³]
- `PsiL50_H`: Leaf water potential at 50% conductance loss for high vegetation [MPa]
- `PsiL50_L`: Leaf water potential at 50% conductance loss for low vegetation [MPa]
- `PsiX50_H`: Xylem water potential at 50% conductance loss for high vegetation [MPa]
- `PsiX50_L`: Xylem water potential at 50% conductance loss for low vegetation [MPa]

# Returns
- `V`: Updated soil water volume [mm]
- `O`: Volumetric soil moisture content [m³/m³]
- `OS`: Surface soil moisture [m³/m³]
- `Psi_soil`: Soil water potential [MPa]
- `Psi_s_H`: Soil water potential for high vegetation [MPa]
- `Psi_s_L`: Soil water potential for low vegetation [MPa]
- `Exwat_H`: Water extraction rate for high vegetation [mm/h]
- `Exwat_L`: Water extraction rate for low vegetation [mm/h]
- `Ko`: Hydraulic conductivity [mm/h]

"""
function soil_moisture_conductivity_update(
    V::Vector{FT},
    Pcla::FT,
    Psan::FT,
    Porg::FT,
    Kfc::FT,
    Phy::FT,
    SPAR::Int,
    Kbot::FT,
    CASE_ROOT_H::Int,
    CASE_ROOT_L::Int,
    ZR95_H::Vector{FT},
    ZR95_L::Vector{FT},
    ZR50_H::Vector{FT},
    ZR50_L::Vector{FT},
    ZRmax_H::Vector{FT},
    ZRmax_L::Vector{FT},
    Zs::Vector{FT},
    Rrootl_H::Vector{FT},
    Rrootl_L::Vector{FT},
    PsiL50_H::Vector{FT},
    PsiL50_L::Vector{FT},
    PsiX50_H::Vector{FT},
    PsiX50_L::Vector{FT},
) where {FT<:AbstractFloat}

    # Calculate soil parameters depending on soil composition
    Zs, dz, ms, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, EvL_Zs, Inf_Zs, RfH_Zs, RfL_Zs, _, Kbot = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT_H,
        CASE_ROOT_L,
        ZR95_H,
        ZR95_L,
        ZR50_H,
        ZR50_L,
        ZRmax_H,
        ZRmax_L,
        Zs,
    )

    # Update soil moisture content in different soil layers
    O, _, _, OS, Psi_s_H, Psi_s_L, _, _, Exwat_H, Exwat_L = soil_water_multilayer(
        V,
        Zs,
        dz,
        ms,
        Osat,
        Ohy,
        nVG,
        alpVG,
        Ks_Zs,
        L,
        Pe,
        O33,
        SPAR,
        EvL_Zs,
        Inf_Zs,
        RfH_Zs,
        RfL_Zs,
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    # Calculate hydraulic conductivity and soil water potential
    Ko, Psi_soil = conductivity_suction(SPAR, Ks_Zs, Osat, Ohy, L, Pe, O33, alpVG, nVG, O)

    return V, O, OS, Psi_soil, Psi_s_H, Psi_s_L, Exwat_H, Exwat_L, Ko
end
