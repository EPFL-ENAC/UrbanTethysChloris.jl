"""
    water_vegetation(
        Rain::FT,
        E_veg::FT,
        In_veg_tm1::FT,
        Sp_In::FT,
        LAI::FT,
        SAI::FT,
        row::FT,
        dth::FT
    ) where {FT<:AbstractFloat}

Calculate canopy rain interception and water balance for vegetation.

# Arguments
- `Rain`: precipitation [mm/dth]
- `E_veg`: evaporation from vegetation [kg/m²s]
- `In_veg_tm1`: Interception from previous time step [mm]
- `Sp_In`: Specific water interception capacity [mm]
- `LAI`: Leaf area index [-]
- `SAI`: Stem area index [-]
- `row`: density of water [kg/m³]
- `dth`: calculation time step [h]

# Returns
- `q_runon_veg`: Runoff [mm/dth]
- `In_veg`: Interception [mm]
- `dIn_veg_dt`: Change in interception [mm/dth]
- `WBalance_In_veg`: Water balance check [mm/dth]
"""
function water_vegetation(
    Rain::FT, E_veg::FT, In_veg_tm1::FT, Sp_In::FT, LAI::FT, SAI::FT, row::FT, dth::FT
) where {FT<:AbstractFloat}
    # Canopy parameters
    PAI = LAI + SAI
    Kthroughfall = FT(0.75)                          # Ramirez and Senarath(2000)
    Cfol = 1 - exp(-Kthroughfall * PAI)             # [-]
    gc = FT(3.7)                                     # [1/mm]
    Kc = FT(0.001) * 60 * dth                       # [mm/dth] -- Mahfouf and Jacquemin 1989
    In_max_veg = Sp_In * (LAI + SAI)                # Sp = [mm]

    # One vegetation layer
    Rain_fol = Rain * Cfol                          # [mm/dth]
    Rain_tf = Rain * (1 - Cfol)                     # [mm/dth]

    # Water balance calculations
    In_veg = In_veg_tm1 + Rain_fol - E_veg * dth * 3600 * 1000 / row  # [mm]
    SE_veg = (In_veg > In_max_veg) * (In_veg - In_max_veg)            # [mm] Storage Excess
    In_veg = In_veg - SE_veg                                           # [mm]

    Dr_veg = Kc * exp(gc * (In_veg - In_max_veg)) * (In_veg > 0)      # [mm/dth] Drainage first layer
    In_veg = In_veg - Dr_veg                                           # [mm]
    Dr_veg = Dr_veg + In_veg * (In_veg < 0)                           # [mm/dth]
    In_veg = max(In_veg, zero(FT))                                     # [mm] First updated Interception

    q_runon_veg = Dr_veg + SE_veg + Rain_tf                           # [mm/dth] Dripping, Saturation Excess and Throughfall
    dIn_veg_dt = In_veg - In_veg_tm1                                  # [mm/dth]

    # Volume Balance check
    WBalance_In_veg = Rain - E_veg * dth * 3600 * 1000 / row - q_runon_veg - dIn_veg_dt  # [mm/dth]

    return q_runon_veg, In_veg, dIn_veg_dt, WBalance_In_veg
end
