"""
    soil_water_multilayer!(hav, hsv, vav, ps, af, timestep, day)

Calculate soil water dynamics and potential water extraction rates.

# Arguments
- `hav`: Hydrologic auxiliary variables
- `hsv`: Hydrologic state variables
- `vav`: Vegetation auxiliary variables
- `ps`: Model parameter set
- `af`: Anthropogenic forcing inputs
- `timestep`: Current simulation timestep
- `day`: Current simulation day

Returns a tuple (Exwat_H, Exwat_L) where:
- Exwat_H: Matrix [cc × n] of potential water extraction rates for high vegetation [mm/h]
- Exwat_L: Matrix [cc × n] of potential water extraction rates for low vegetation [mm/h]
where cc is the number of crown areas and n is the number of soil layers.
"""
function soil_water_multilayer(
    V::Vector{FT},
    Zs::Vector{FT},
    dz::Vector{FT},
    n::Int,
    Osat::Vector{FT},
    Ohy::Vector{FT},
    nVG::Vector{FT},
    alpVG::Vector{FT},
    Ks_Zs::Vector{FT},
    L::Vector{FT},
    Pe::Vector{FT},
    O33::Vector{FT},
    SPAR::Int,
    EvL_Zs::Vector{FT},
    Inf_Zs::Vector{FT},
    RfH_Zs::Matrix{FT},
    RfL_Zs::Matrix{FT},
    Rrootl_H::Vector{FT},
    Rrootl_L::Vector{FT},
    PsiL50_H::Vector{FT},
    PsiL50_L::Vector{FT},
    PsiX50_H::Vector{FT},
    PsiX50_L::Vector{FT},
) where {FT<:AbstractFloat}
    O = ones(FT, n)
    WTR = zeros(FT, n)
    for i in n:-1:1
        if i == n
            O[i] = (V[i] / dz[i]) + Ohy[i]  # Water Content [] All layers
            WTR[i] = (O[i] - Osat[i]) * dz[i] * (O[i] > Osat[i])  # [mm] Water Table Rise
        else
            O[i] = (V[i] + WTR[i + 1]) / dz[i] + Ohy[i]  # Water Content [] All layers
            WTR[i] = (O[i] - Osat[i]) * dz[i] * (O[i] > Osat[i])  # [mm] Excess From Reservoir - Below
        end
    end

    O = clamp.(O, Ohy, Osat)

    # Find water table depth
    i = n
    while i > 0 && O[i] > (Osat[i] - 1e-5)
        i -= 1
    end

    ZWT = Zs[i + 1]
    PHead = @. (Zs[1:n] + dz / 2) - ZWT
    PHead = max.(PHead, 0)

    # Calculate potential based on SPAR
    POT = compute_potential(Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, O, PHead)

    # Recalculate water table depth
    if i == n
        ZWT = Zs[i + 1]
    elseif i == 0
        ZWT = 0
    else
        if POT[i + 1] > 0
            CZ = Zs[1:n] .+ dz / 2
            ZWT = CZ[i + 1] + (CZ[i] - CZ[i + 1]) * (-POT[i + 1]) / (POT[i] - POT[i + 1])
        else
            ZWT = Zs[i + 1]
        end

        PHead = @. (Zs[1:n] + dz / 2) - ZWT
        PHead = max.(PHead, 0)

        POT = compute_potential(Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, O, PHead)
    end

    Rd = WTR[1]  # [mm] Dunne Runoff
    OS = sum(EvL_Zs .* O)  # Evaporation Bare Soil WC []
    OF = sum(Inf_Zs .* O)  # Infiltration Water Content []

    cc = size(RfH_Zs, 1)
    OH = zeros(FT, cc)
    OL = zeros(FT, cc)
    Psi_s_H = zeros(FT, cc)
    Psi_s_L = zeros(FT, cc)

    for i in 1:cc
        OH[i] = sum(view(RfH_Zs, i, :) .* O)
        OL[i] = sum(view(RfL_Zs, i, :) .* O)

        _, Psi_s_H[i] = weighted_conductivity_suction(
            SPAR, Ks_Zs, Osat, Ohy, L, Pe, O33, alpVG, nVG, OH[i], RfH_Zs[i, :]
        )

        _, Psi_s_L[i] = weighted_conductivity_suction(
            SPAR, Ks_Zs, Osat, Ohy, L, Pe, O33, alpVG, nVG, OL[i], RfL_Zs[i, :]
        )
    end

    # Convert to MPa
    Psi_s_H = @. -(Psi_s_H / 1000) * 1000 * 9.81 / 1e6
    Psi_s_L = @. -(Psi_s_L / 1000) * 1000 * 9.81 / 1e6

    rho2 = 55555555.0  # [mmolH20 /m^3] Water density
    rcyl = 2e-3        # [m] radius cylinder of soil to which root has access to
    rroot = 0.5e-3     # [m] root radius

    Psi_s = zeros(FT, n)
    gsr_L = zeros(FT, cc, n)
    gsr_H = zeros(FT, cc, n)

    # Layer by layer analysis

    Ko, Psi_s = conductivity_suction(SPAR, Ks_Zs, Osat, Ohy, L, Pe, O33, alpVG, nVG, O)

    for jk in 1:n
        for i in 1:cc
            gsr_L[i, jk] = root_soil_conductance(
                Ko[jk], RfL_Zs[i, jk] * Rrootl_L[i], rcyl, rroot, dz[jk]
            )
            gsr_H[i, jk] = root_soil_conductance(
                Ko[jk], RfH_Zs[i, jk] * Rrootl_H[i], rcyl, rroot, dz[jk]
            )
        end
    end

    Psi_s = @. -(Psi_s / 1000) * 1000 * 9.81 / 1e6

    Psi_minH = minimum([PsiX50_H; PsiL50_H]; dims=1)
    Psi_minL = minimum([PsiX50_L; PsiL50_L]; dims=1)

    Exwat_L = gsr_L / rho2 * 1000.0 * 3600.0 .* (-Psi_minL' .+ Psi_s')
    Exwat_H = gsr_H / rho2 * 1000.0 * 3600.0 .* (-Psi_minH' .+ Psi_s')

    Exwat_L[Exwat_L .< 0] .= 0.0
    Exwat_H[Exwat_H .< 0] .= 0.0

    gsr_L = sum(gsr_L; dims=2)
    gsr_H = sum(gsr_H; dims=2)

    return O,
    ZWT, OF, OS, Psi_s_H, Psi_s_L, gsr_H, gsr_L, Exwat_H, Exwat_L, Rd, WTR, POT, OH,
    OL
end

function compute_potential(
    Osat::Vector{FT},
    Ohy::Vector{FT},
    nVG::Vector{FT},
    alpVG::Vector{FT},
    Ks_Zs::Vector{FT},
    L::Vector{FT},
    Pe::Vector{FT},
    O33::Vector{FT},
    SPAR::Int,
    O::Vector{FT},
    PHead::Vector{FT},
) where {FT<:AbstractFloat}
    POT = if SPAR == 1
        Se = @. (O - Ohy) / (Osat - Ohy)
        mVG = @. 1 - 1 / nVG
        PHead + (1.0 ./ alpVG) * ((Se) .^ (-1.0 ./ mVG) - 1.0) .^ (1.0 ./ nVG)
    elseif SPAR == 2
        _, Ptem = conductivity_suction(SPAR, Ks_Zs, Osat, Ohy, L, Pe, O33, alpVG, nVG, O)
        PHead - Ptem
    end

    return POT
end
