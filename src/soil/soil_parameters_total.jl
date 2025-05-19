"""
    soil_parameters_total(
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
        Zs::Vector{FT}
    ) where {FT<:AbstractFloat}

Calculate total soil parameters based on composition and root distribution.

# Arguments
- `Pcla`: Clay fraction in soil [-]
- `Psan`: Sand fraction in soil [-]
- `Porg`: Organic matter fraction in soil [-]
- `Kfc`: Hydraulic conductivity at field capacity [mm/h]
- `Phy`: Soil water potential at hygroscopic point [kPa]
- `SPAR`: Soil parameterization choice:
    1. van Genuchten (1980)
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

# Returns
- `Zs`: Soil layer depths [mm]
- `dz`: Layer thicknesses [mm]
- `ms`: Number of soil layers [-]
- `Osat`: Saturated water content [m³/m³]
- `Ohy`: Hygroscopic water content [m³/m³]
- `nVG`: van Genuchten n parameter [-]
- `alpVG`: van Genuchten α parameter [1/mm]
- `Ks_Zs`: Saturated hydraulic conductivity [mm/h]
- `L`: Pore size distribution index [-]
- `Pe`: Air entry pressure [kPa]
- `O33`: Water content at -33 kPa [m³/m³]
- `SPAR`: Soil parameterization type [-]
- `EvL_Zs`: Evaporation layer fractions [-]
- `Inf_Zs`: Infiltration layer fractions [-]
- `RfH_Zs`: Root fraction distribution for high vegetation [-]
- `RfL_Zs`: Root fraction distribution for low vegetation [-]
- `Zinf`: Infiltration depth [mm]
- `Kbot`: Bottom boundary conductivity [mm/h]
- `Slo_pot`: Slope fractions [-]
- `Dz`: Layer center distances [mm]
- `aR`: Horizontal length scale [-]
- `aTop`: Area to contour length ratio [mm]
- `rsd`: Dry soil density [kg/m³]
- `lan_dry`: Dry soil thermal conductivity [W/m K]
- `lan_s`: Solid soil thermal conductivity [W/m K]
- `cv_s`: Solid soil volumetric heat capacity [J/m³ K]
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
