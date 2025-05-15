"""
    soil_parameters2(
        Osat::Vector{FT},
        L::Vector{FT},
        Pe::Vector{FT},
        Ks::Vector{FT},
        O33::Vector{FT},
        nVG::Vector{FT},
        alpVG::Vector{FT},
        Kfc::FT,
        Pss::FT,
        Pwp::FT,
        Phy::FT,
        SPAR::Int = 2
    ) where {FT<:AbstractFloat}

Calculate soil moisture retention points based on soil hydraulic parameters.

# Arguments
- `Osat::Vector{FT}`: Saturated soil water content
- `L::Vector{FT}`: Pore size distribution index
- `Pe::Vector{FT}`: Air entry pressure
- `Ks::Vector{FT}`: Saturated hydraulic conductivity
- `O33::Vector{FT}`: Soil water content at -33 kPa
- `nVG::Vector{FT}`: van Genuchten n parameter
- `alpVG::Vector{FT}`: van Genuchten α parameter
- `Kfc::FT`: Hydraulic conductivity at field capacity
- `Pss::FT`: Soil water potential at steady state
- `Pwp::FT`: Permanent wilting point potential
- `Phy::FT`: Hygroscopic point potential
- `SPAR::Int=2`: Soil parameter set (1: van Genuchten, 2: Brooks-Corey)

# Returns
A tuple of four vectors:
- `Ofc`: Soil water content at field capacity
- `Oss`: Soil water content at steady state
- `Owp`: Soil water content at wilting point
- `Ohy`: Soil water content at hygroscopic point
"""
function soil_parameters2(
    Osat::Vector{FT},
    L::Vector{FT},
    Pe::Vector{FT},
    Ks::Vector{FT},
    O33::Vector{FT},
    nVG::Vector{FT},
    alpVG::Vector{FT},
    Kfc::FT,
    Pss::FT,
    Pwp::FT,
    Phy::FT,
    Ohy::Vector{FT}=zeros(FT, length(Osat)),
    SPAR::Int=2,
) where {FT<:AbstractFloat}

    # Initialize output arrays
    ms = length(Osat)
    Ofc = zeros(FT, ms)
    Oss = zeros(FT, ms)
    Owp = zeros(FT, ms)

    if SPAR == 1
        Ohy = copy(Ohy)
        mVG = 1 .- 1 ./ nVG
        Pss = -101.9368 * Pss  # [mm]
        Se = 1 ./ ((1 .+ abs.(alpVG * Pss) .^ nVG) .^ mVG)
        Se[Se .> 1] .= 1
        O = Ohy .+ (Osat .- Ohy) .* Se
        Oss[:] = O

        Pwp = -101.9368 * Pwp  # [mm]
        Se = 1 ./ ((1 .+ abs.(alpVG * Pwp) .^ nVG) .^ mVG)
        Se[Se .> 1] .= 1
        O = Ohy .+ (Osat .- Ohy) .* Se
        Owp[:] = O
        Ohy = zeros(FT, ms)

    elseif SPAR == 2
        for i in 1:ms
            B = 1/L[i]
            A = exp(log(33) + B * log(O33[i]))  # Coefficient of moisture tension

            if Pss < 33
                Oss[i] = O33[i] + (33 - Pss) * (Osat[i] - O33[i]) / (33 - Pe[i])
            else
                Oss[i] = (Pss / A)^(-1 / B)
            end

            if Pwp < 33
                Owp[i] = O33[i] + (33 - Pwp) * (Osat[i] - O33[i]) / (33 - Pe[i])
            else
                Owp[i] = (Pwp / A)^(-1 / B)
            end

            if Phy < 33
                Ohy[i] = O33[i] + (33 - Phy) * (Osat[i] - O33[i]) / (33 - Pe[i])
            else
                Ohy[i] = (Phy / A)^(-1 / B)
            end

            Ofc[i] = Osat[i] * (Kfc / Ks[i])^(1 / (3 + 2/L[i]))
        end
    end

    return Ofc, Oss, Owp, Ohy
end
