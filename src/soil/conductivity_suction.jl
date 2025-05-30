"""
    conductivity_suction(
        SPAR::Int,
        Ks::FT,
        Osat::FT,
        Ohy::FT,
        L::FT,
        Pe::FT,
        O33::FT,
        alpVG::FT,
        nVG::FT,
        O::FT
    ) where {FT<:AbstractFloat}

    conductivity_suction(
        SPAR::Int,
        Ks::FT,
        Osat::FT,
        Ohy::FT,
        L::FT,
        Pe::FT,
        O33::FT,
        alpVG::FT,
        nVG::FT,
        O::Vector{FT},
    ) where {FT<:AbstractFloat}

    conductivity_suction(
        SPAR::Int,
        Ks::Vector{FT},
        Osat::Vector{FT},
        Ohy::Vector{FT},
        L::Vector{FT},
        Pe::Vector{FT},
        O33::Vector{FT},
        alpVG::Vector{FT},
        nVG::Vector{FT},
        O::Vector{FT},
    ) where {FT<:AbstractFloat}

    weighted_conductivity_suction(
        SPAR::Int,
        Ks::Vector{FT},
        Osat::Vector{FT},
        Ohy::Vector{FT},
        L::Vector{FT},
        Pe::Vector{FT},
        O33::Vector{FT},
        alpVG::Vector{FT},
        nVG::Vector{FT},
        O::FT,
        weights::Vector{FT},
    ) where {FT<:AbstractFloat}

Calculate soil hydraulic conductivity and matric suction using specified model.

# Arguments
- `SPAR`: Soil parameterization choice:
    1. Van-Genuchten (1980) corrected
    2. Saxton-Rawls (1986)
- `Ks`: Saturated hydraulic conductivity [mm/h]
- `Osat`: Saturated water content [m³/m³]
- `Ohy`: Hygroscopic water content [m³/m³]
- `L`: Lambda parameter for Saxton-Rawls model [-]
- `Pe`: Air entry pressure [kPa]
- `O33`: Water content at -33 kPa [m³/m³]
- `alpVG`: Van Genuchten α parameter [1/mm]
- `nVG`: Van Genuchten n parameter [-]
- `O`: Current water content [m³/m³]
- `weights`: Weights for each soil layer [dimensionless]

# Returns
- `Ko`: Unsaturated hydraulic conductivity [mm/h]
- `Po`: Matric potential [mm]
"""
function conductivity_suction(
    SPAR::Int, Ks::FT, Osat::FT, Ohy::FT, L::FT, Pe::FT, O33::FT, alpVG::FT, nVG::FT, O::FT
) where {FT<:AbstractFloat}
    if SPAR == 1
        # Van-Genuchten, 1980 Corrected
        Se = (O .- Ohy) ./ (Osat .- Ohy)
        mVG = 1 .- 1 ./ nVG
        Po = -(1 ./ alpVG) .* ((Se) .^ (-1 ./ mVG)-1) .^ (1 ./ nVG)
        Ko = Ks .* ((Se) .^ (0.5)) .* (1 .- (1 .- Se .^ (1 ./ mVG)) .^ mVG) .^ 2  # [mm/h]

    elseif SPAR == 2
        # Saxton and Rawls 1986
        gw = 9810.0  # specific weight water [N/m^3]
        B = 1 ./ L
        A = exp.(log(33) .+ B .* log.(O33))  # Coefficient of moisture tension

        Ko = Ks .* (O ./ Osat) .^ (3 .+ 2 ./ L)  # [mm/h]

        if O < O33
            Psi = A .* O .^ (-B)  # [kPa]
        else
            Psi = 33 .- ((O .- O33) .* (33 .- Pe) ./ (Osat .- O33))  # [kPa]
        end
        Po = 1000 * 1000 * Psi ./ gw  # [mm]
    end

    return Ko, Po
end

function conductivity_suction(
    SPAR::Int,
    Ks::FT,
    Osat::FT,
    Ohy::FT,
    L::FT,
    Pe::FT,
    O33::FT,
    alpVG::FT,
    nVG::FT,
    O::Vector{FT},
) where {FT<:AbstractFloat}
    Ko = zeros(FT, length(O))
    Po = zeros(FT, length(O))

    for i in 1:length(O)
        Ko[i], Po[i] = conductivity_suction(
            SPAR, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, O[i]
        )
    end

    return Ko, Po
end

function conductivity_suction(
    SPAR::Int,
    Ks::Vector{FT},
    Osat::Vector{FT},
    Ohy::Vector{FT},
    L::Vector{FT},
    Pe::Vector{FT},
    O33::Vector{FT},
    alpVG::Vector{FT},
    nVG::Vector{FT},
    O::Vector{FT},
) where {FT<:AbstractFloat}
    Ko = zeros(FT, length(O))
    Po = zeros(FT, length(O))

    for i in 1:length(O)
        Ko[i], Po[i] = conductivity_suction(
            SPAR, Ks[i], Osat[i], Ohy[i], L[i], Pe[i], O33[i], alpVG[i], nVG[i], O[i]
        )
    end

    return Ko, Po
end

function weighted_conductivity_suction(
    SPAR::Int,
    Ks::Vector{FT},
    Osat::Vector{FT},
    Ohy::Vector{FT},
    L::Vector{FT},
    Pe::Vector{FT},
    O33::Vector{FT},
    alpVG::Vector{FT},
    nVG::Vector{FT},
    O::FT,
    weights::Vector{FT},
) where {FT<:AbstractFloat}
    Ks_weighted = sum(Ks .* weights)
    Osat_weighted = sum(Osat .* weights)
    Ohy_weighted = sum(Ohy .* weights)
    L_weighted = sum(L .* weights)
    Pe_weighted = sum(Pe .* weights)
    O33_weighted = sum(O33 .* weights)
    alpVG_weighted = sum(alpVG .* weights)
    nVG_weighted = sum(nVG .* weights)

    return conductivity_suction(
        SPAR,
        Ks_weighted,
        Osat_weighted,
        Ohy_weighted,
        L_weighted,
        Pe_weighted,
        O33_weighted,
        alpVG_weighted,
        nVG_weighted,
        O,
    )
end
