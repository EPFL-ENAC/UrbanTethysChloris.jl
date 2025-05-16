"""
    leakage_bottom(O, Ks_Zs, Osat, Ohy, L, nVG, Kbot, ms, SPAR)

Calculate bottom leakage for soil water movement.

# Arguments
- `O`: Soil water content
- `Ks_Zs`: Saturated hydraulic conductivity [mm/h]
- `Osat`: Saturated water content
- `Ohy`: Residual water content
- `L`: Lambda parameter
- `nVG`: van Genuchten n parameter
- `Kbot`: Bottom hydraulic conductivity [mm/h]
- `ms`: Layer index
- `SPAR`: Soil parameter set selector

# Returns
- `Lk`: Bottom leakage [mm/h]
"""
function leakage_bottom(
    O::Vector{FT},
    Ks_Zs::Vector{FT},
    Osat::Vector{FT},
    Ohy::Vector{FT},
    L::Vector{FT},
    nVG::Vector{FT},
    Kbot::FT,
    ms::Int,
    SPAR::Int,
) where {FT<:AbstractFloat}
    if isnan(Kbot)
        if SPAR == 1
            Se = (O[ms] - Ohy[ms]) / (Osat[ms] - Ohy[ms])
            mVG = 1 - 1/nVG[ms]
            Ko = Ks_Zs[ms] * (Se^0.5) * (1 - (1 - Se^(1/mVG))^mVG)^2
        elseif SPAR == 2
            Ko = Ks_Zs[ms] * (O[ms]/Osat[ms])^(3 + 2/L[ms])
        end

        return Ko
    else
        if O[ms] > Osat[ms] - 1e-5
            return Kbot
        else
            return zero(FT)
        end
    end
end
