"""
    soil_moisture_conductivity_update(
        V::Vector{FT}, Pcla::FT, Psan::FT, Porg::FT, Kfc::FT, Phy::FT, SPAR::Int,
        Kbot::FT, CASE_ROOT_H::Int, CASE_ROOT_L::Int, ZR95_H::Vector{FT}, ZR95_L::Vector{FT},
        ZR50_H::Vector{FT}, ZR50_L::Vector{FT}, ZRmax_H::Vector{FT}, ZRmax_L::Vector{FT},
        Zs::Vector{FT}, Rrootl_H::Vector{FT}, Rrootl_L::Vector{FT}, PsiL50_H::Vector{FT},
        PsiL50_L::Vector{FT}, PsiX50_H::Vector{FT}, PsiX50_L::Vector{FT}
    ) where {FT<:AbstractFloat}

Updates soil moisture content and calculates hydraulic conductivity and soil water potential.

Returns a tuple containing:
- `V`: Updated soil water volume
- `O`: Soil moisture content
- `OS`: Surface soil moisture
- `Psi_soil`: Soil water potential
- `Psi_s_H`: Soil water potential for high vegetation
- `Psi_s_L`: Soil water potential for low vegetation
- `Exwat_H`: Water extraction rate for high vegetation
- `Exwat_L`: Water extraction rate for low vegetation
- `Ko`: Hydraulic conductivity
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
