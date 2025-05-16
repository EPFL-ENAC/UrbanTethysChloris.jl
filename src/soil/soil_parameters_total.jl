"""
    soil_parameters_total(
        Pcla::FT,      # Fraction of clay in the soil [-]
        Psan::FT,      # Fraction of sand in the soil [-]
        Porg::FT,      # Fraction of organic material in the soil [-]
        Kfc::FT,       # Conductivity at field capacity [mm/h]
        Phy::FT,       # Suction at the residual/hygroscopic water content [kPa]
        SPAR::Int,     # SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
        Kbot::FT,      # Conductivity at the bedrock layer [mm/h]
        CASE_ROOT_H::Int,  # Type of Root Profile for high vegetation
        CASE_ROOT_L::Int,  # Type of Root Profile for low vegetation
        ZR95_H::Vector{FT},  # Root depth 95 percentile for high vegetation [mm]
        ZR95_L::Vector{FT},  # Root depth 95 percentile for low vegetation [mm]
        ZR50_H::Vector{FT},  # Root depth 50 percentile for high vegetation [mm]
        ZR50_L::Vector{FT},  # Root depth 50 percentile for low vegetation [mm]
        ZRmax_H::Vector{FT}, # Maximum Root depth for high vegetation [mm]
        ZRmax_L::Vector{FT}, # Maximum Root depth for low vegetation [mm]
        Zs::Vector{FT}       # soil layer discretization [mm]
    ) where {FT<:AbstractFloat}

Returns:
- `Zs::Vector{FT}`: soil layer discretization [mm]
- `dz::Vector{FT}`: Thickness of soil layers [mm]
- `ms::Int`: Number of soil layers [-]
- `Osat::Vector{FT}`: Saturation moisture 0 kPa [-]
- `Ohy::Vector{FT}`: Hygroscopic Moisture at 10000 kPa [-]
- `nVG::Vector{FT}`: n parameter Van-Genuchten curve [1/mm]
- `alpVG::Vector{FT}`: Alpha parameter Van-Genuchten curve [1/mm]
- `Ks_Zs::Vector{FT}`: Hydraulic conductivity at saturation [mm/h]
- `L::Vector{FT}`: Slope of logaritimc tension-moisture curve [-]
- `Pe::Vector{FT}`: Tension at air entry (bubbling pressure) [kPa]
- `O33::Vector{FT}`: 33 kPa Moisture [-]
- `SPAR::Int`: SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
- `EvL_Zs::Vector{FT}`: Fraction of evaporation depth per layer [-]
- `Inf_Zs::Vector{FT}`: Fraction of infiltration depth per layer [-]
- `RfH_Zs::Matrix{FT}`: Root Fraction for high vegetation [-]
- `RfL_Zs::Matrix{FT}`: Root Fraction for low vegetation [-]
- `Zinf::FT`: Depth of infiltration layer [mm]
- `Kbot::FT`: Conductivity at bedrock layer [mm/h]
- `Slo_pot::Vector{FT}`: Fraction dy/dx [-]
- `Dz::Vector{FT}`: Delta depth between layers [mm]
- `aR::FT`: anisotropy ratio [-]
- `aTop::FT`: Ratio Area/ContourLength [mm]
- `rsd::Vector{FT}`: Normal density dry soil [kg/m^3]
- `lan_dry::Vector{FT}`: Thermal conductivity dry soil [W/m K]
- `lan_s::Vector{FT}`: Thermal conductivity solid soil [W/m K]
- `cv_s::Vector{FT}`: Volumetric heat capacity solid soil [J/m^3 K]
"""
function soil_parameters_total(
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
) where {FT<:AbstractFloat}
    ms = length(Zs) - 1
    dz = diff(Zs)
    Dz = zeros(FT, ms)

    for i in 1:ms
        if i > 1
            Dz[i] = (dz[i] + dz[i - 1])/2  # Delta Depth Between Middle Layer [mm]
        else
            Dz[i] = dz[1]/2  # Delta Depth Between First Middle Layer and soil surface [mm]
        end
    end

    Zdes = Zs[2] - Zs[1]  # Depth of desorption/evaporation layer
    Zinf = Zs[2] - Zs[1]  # Depth of infiltration layer

    EvL_Zs = evaporation_layers(Zs, Zdes)  # Fraction of evaporation depth
    Inf_Zs = evaporation_layers(Zs, Zinf)  # Fraction of infiltration depth

    Slo_pot = zeros(FT, ms)  # [fraction dy/dx]
    aR = one(FT)  # anisotropy ratio
    cellsize = one(FT)
    aTop = 1000*cellsize^2/cellsize  # [mm] Ratio between Area/ContourLength

    # Soil Parameters from composition
    Osat, L, Pe, Ks, O33, rsd, lan_dry, lan_s, cv_s = soil_parameters(Psan, Pcla, Porg)

    # Convert scalar parameters to vectors
    rsd = fill(rsd, ms)
    lan_dry = fill(lan_dry, ms)
    lan_s = fill(lan_s, ms)
    cv_s = fill(cv_s, ms)

    # Van Genuchten parameters
    p = 3 + 2/L
    m = 2/(p-1)
    nVG = 1/(1-m)
    alpVG =
        (((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1

    # Initialize vectors for each soil layer
    Osat = fill(Osat, ms)
    L = fill(L, ms)
    Pe = fill(Pe, ms)
    Ks_Zs = fill(Ks, ms)
    O33 = fill(O33, ms)
    nVG = fill(nVG, ms)
    alpVG = fill(alpVG, ms)

    # Soil parameters II
    _, _, _, Ohy = soil_parameters2(Osat, L, Pe, Ks_Zs, O33, nVG, alpVG, Kfc, NaN, NaN, Phy)

    # Root Distribution
    if CASE_ROOT_H != CASE_ROOT_L
        @warn "CASE_ROOT_H and CASE_ROOT_L are not the same. Using CASE_ROOT_H."
    end

    RfH_Zs, RfL_Zs = root_fraction_general(
        Zs, CASE_ROOT_H, ZR95_H, ZR50_H, ZR95_L, ZR50_L, ZRmax_H, ZRmax_L
    )

    return Zs,
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
    Zinf,
    Kbot,
    Slo_pot,
    Dz,
    aR,
    aTop,
    rsd,
    lan_dry,
    lan_s,
    cv_s
end
