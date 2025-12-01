"""
    soil_moistures_rich_comp(
        V,
        Lk,
        f,
        EG,
        T_H,
        T_L,
        Qi_in,
        IS,
        SPAR,
        Osat,
        Ohy,
        O33,
        dz,
        Ks_Zs,
        Dz,
        numn,
        L,
        Pe,
        aR,
        aT,
        alpVG,
        nVG,
        cosalp,
        sinalp,
        SN)

Calculate soil moisture changes using Richards equation-based hydrological model.

# Arguments
- `V::Vector{FT}`: Current soil water volume per unit area [mm]
- `Lk::FT`: Bottom leakage rate [mm/h]
- `f::FT`: Infiltration rate at top boundary [mm/h]
- `EG::Vector{FT}`: Evaporation rates from ground [mm/h]
- `T_H::Vector{FT}`: Transpiration rates from hydrophilic roots [mm/h]
- `T_L::Vector{FT}`: Transpiration rates from hydrophobic roots [mm/h]
- `Qi_in::Vector{FT}`: Lateral inflow rates [mm/h]
- `IS::Vector{FT}`: Surface water indicator (1 for surface water present)
- `SPAR::Int`: Soil parameter set selector (1: Van Genuchten, 2: Brooks-Corey)
- `Osat::Vector{FT}`: Saturated water content at 0 kPa
- `Ohy::Vector{FT}`: Residual (hygroscopic) water content
- `O33::Vector{FT}`: Water content at -33 kPa tension
- `dz::Vector{FT}`: Layer thicknesses [mm]
- `Ks_Zs::Vector{FT}`: Saturated hydraulic conductivity [mm/h]
- `Dz::Vector{FT}`: Distance between layer centers [mm]
- `numn::Int`: Number of soil layers
- `L::Vector{FT}`: Lambda parameter (slope of log tension-moisture curve)
- `Pe::Vector{FT}`: Air entry tension (bubbling pressure) [kPa]
- `aR::FT`: Horizontal length scale parameter
- `aT::FT`: Terrain curvature parameter
- `alpVG::Vector{FT}`: van Genuchten α parameter
- `nVG::Vector{FT}`: van Genuchten n parameter
- `cosalp::FT`: Cosine of slope angle
- `sinalp::FT`: Sine of slope angle
- `SN::FT`: Stream network identifier

# Returns
- `dV::Vector{FT}`: Change rates of soil water volume per unit area [mm/h]
"""
function soil_moistures_rich_comp(
    V::Vector{FT},
    Lk::FT,
    f::FT,
    EG::Vector{FT},
    T_H::Vector{FT},
    T_L::Vector{FT},
    Qi_in::Vector{FT},
    IS::Vector{FT},
    SPAR::Int, # TODO: replace with enum-like type
    Osat::Vector{FT},
    Ohy::Vector{FT},
    O33::Vector{FT},
    dz::Vector{FT},
    Ks_Zs::Vector{FT},
    Dz::Vector{FT},
    numn::Int,
    L::Vector{FT},
    Pe::Vector{FT},
    aR::FT,
    aT::FT,
    alpVG::Vector{FT},
    nVG::Vector{FT},
    cosalp::FT,
    sinalp::FT,
    SN::FT,
) where {FT<:AbstractFloat}
    O = V ./ dz .+ Ohy
    I1 = O .>= Osat .- 1e-5
    I2 = O .<= Ohy .+ 1e-5
    O[I1] .= Osat[I1]
    O[I2] .= Ohy[I2] .+ 1e-5

    numnm1 = numn - 1

    dV = zeros(FT, numn)
    K = zeros(FT, numn)
    P = zeros(FT, numn)
    Khalf = zeros(FT, numnm1)
    q = zeros(FT, numnm1)
    qf = zeros(FT, numnm1)

    if SPAR == 1
        mVG = 1 .- 1 ./ nVG
        Se = (O .- Ohy) ./ (Osat .- Ohy)
        P = (1 ./ alpVG) .* ((Se) .^ (-1 ./ mVG) .- 1) .^ (1 ./ nVG)
        K = Ks_Zs .* ((Se) .^ 0.5) .* (1 .- (1 .- (Se) .^ (1 ./ mVG)) .^ mVG) .^ 2
    elseif SPAR == 2
        B = 1 ./ L
        A = exp.(log.(33) .+ B .* log.(O33))

        for i in 1:numn
            K[i] = Ks_Zs[i] * (O[i] / Osat[i])^(3 + (2 / L[i]))
            if O[i] < O33[i]
                P[i] = A[i] * (O[i]^-B[i])
            else
                P[i] = 33 - ((O[i] - O33[i]) * (33 - Pe[i]) / (Osat[i] - O33[i]))
            end
        end
        P .= -101.9368 .* P
    end

    To = K .* aR .* dz

    Qi_out = (To ./ aT) .* sinalp .* (SN == 1 ? IS : 1.0)

    Khalf[1:numnm1] = 0.5 .* (K[1:numnm1] .+ K[2:numn])
    q = Khalf .* (1 * cosalp .- (P[2:numn] .- P[1:numnm1]) ./ Dz[2:end])
    qf = (1 .- I1[2:numn]) .* q
    q[q .> 0] .= qf[q .> 0]

    dV[1] = f - q[1] - T_H[1] - T_L[1] - EG[1] + Qi_in[1] - Qi_out[1]
    for i in 2:(numn - 1)
        dV[i] = q[i - 1] - q[i] - T_H[i] - T_L[i] - EG[i] + Qi_in[i] - Qi_out[i]
    end
    dV[numn] =
        q[numnm1] - Lk - T_H[numn] - T_L[numn] - EG[numn] + Qi_in[numn] - Qi_out[numn]

    return dV
end
