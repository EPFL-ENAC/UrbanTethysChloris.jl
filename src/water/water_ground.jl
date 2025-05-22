"""
    water_ground(
        q_runon_veg::FT,
        Runon_tm1::FT,
        E_ground::FT,
        Otm1::Vector{FT},
        In_ground_tm1::FT,
        In_max_ground::FT,
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
        dth::FT,
        row::FT
    ) where {FT<:AbstractFloat}

Calculate ground water dynamics and potential infiltration rates.

# Arguments
- `q_runon_veg`: Runon from vegetation [mm/dth]
- `Runon_tm1`: Previous timestep runon [mm/dth]
- `E_ground`: Ground evaporation [kg/m²s]
- `Otm1`: Previous timestep soil moisture [m³/m³]
- `In_ground_tm1`: Previous timestep ground interception [mm]
- `In_max_ground`: Maximum ground interception capacity [mm]
- `Pcla`: Clay fraction in soil [-]
- `Psan`: Sand fraction in soil [-]
- `Porg`: Organic matter fraction in soil [-]
- `Kfc`: Hydraulic conductivity at field capacity [mm/h]
- `Phy`: Soil water potential at hygroscopic point [kPa]
- `SPAR`: Soil parameterization choice [-]
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
- `dth`: Calculation time step [h]
- `row`: Water density [kg/m³]

# Returns
- `q_runon_ground`: Ground runoff [mm/dth]
- `In_ground`: Ground interception [mm]
- `dIn_ground_dt`: Change in ground interception [mm/dth]
- `f_ground`: Infiltration rate [mm/h]
- `WBalance_In_ground`: Water balance check [mm/dth]
"""
function water_ground(
    q_runon_veg::FT,
    Runon_tm1::FT,
    E_ground::FT,
    Otm1::Vector{FT},
    In_ground_tm1::FT,
    In_max_ground::FT,
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
    dth::FT,
    row::FT,
) where {FT<:AbstractFloat}
    # Calculate soil parameters
    _, _, _, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, _, _, _, _, Zinf = Soil.soil_parameters_total(
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

    # Water interception on ground under vegetation
    q_runon_veg = q_runon_veg + Runon_tm1  # [mm/dth]
    In_ground = In_ground_tm1 + q_runon_veg - E_ground * dth * 3600 * 1000 / row  # [mm]
    WIS = In_ground / dth  # [mm/h] total Water Incoming to Soil Layer
    ydepth = In_ground_tm1 / dth  # [mm/h] Interception/ponding from previous time step

    # Calculate infiltration rate [mm/h]
    f_ground, _ = Soil.infiltration(
        Osat[1],
        Ohy[1],
        L[1],
        alpVG[1],
        nVG[1],
        Pe[1],
        Ks_Zs[1],
        O33[1],
        SPAR,
        Otm1[1],
        Zinf,
        WIS,
        one(FT),
        ydepth,
    )

    # Water balance
    In_ground =
        In_ground_tm1 + q_runon_veg - E_ground * dth * 3600 * 1000 / row - f_ground * dth  # [mm]
    q_runon_ground = (In_ground > In_max_ground) * (In_ground - In_max_ground)  # [mm/dth]
    In_ground = In_ground - (In_ground > In_max_ground) * (In_ground - In_max_ground)  # [mm]
    dIn_ground_dt = In_ground - In_ground_tm1  # [mm/dth]

    # Volume Balance check
    WBalance_In_ground =
        q_runon_veg - E_ground * dth * 3600 * 1000 / row - f_ground * dth - q_runon_ground -
        dIn_ground_dt  # [mm/dth]

    return q_runon_ground, In_ground, dIn_ground_dt, f_ground, WBalance_In_ground
end
