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

# Returns
Tuple containing:
- `Ko`: Unsaturated hydraulic conductivity [mm/h]
- `Po`: Matric potential [mm]

# References
- van Genuchten (1980): Soil Science Society of America Journal, 44(5), 892-898
- Saxton & Rawls (1986): Soil Science Society of America Journal, 50(4), 1031-1036
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
