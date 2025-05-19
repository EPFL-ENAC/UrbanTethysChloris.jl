"""
    water_impervious(
        Rain::FT,
        Runon_tm1::FT,
        E_imp::FT,
        In_imp_tm1::FT,
        dth::FT,
        row::FT,
        In_max_imp::FT,
        K_imp::FT
    ) where {FT<:AbstractFloat}

Calculates water balance for impervious surfaces.

# Arguments
- `Rain`: precipitation [mm/dth]
- `Runon_tm1`: Previous timestep runon [mm/dth]
- `E_imp`: evaporation from impervious surfaces [kg/m²s]
- `In_imp_tm1`: Interception from previous time step [mm]
- `dth`: calculation time step [h]
- `row`: density of water [kg/m³]
- `In_max_imp`: Maximum interception capacity of urban area [mm]
- `K_imp`: Hydraulic conductivity of impervious area [mm/h]

# Returns
- `q_runon_imp`: Runoff [mm/dth]
- `In_imp`: Interception [mm]
- `dIn_imp_dt`: Change in interception [mm/dth]
- `Lk_imp`: Leakage from impervious area [mm/h]
- `WBalance_In_imp`: Water balance check [mm/dth]
"""
function water_impervious(
    Rain::FT,
    Runon_tm1::FT,
    E_imp::FT,
    In_imp_tm1::FT,
    dth::FT,
    row::FT,
    In_max_imp::FT,
    K_imp::FT,
) where {FT<:AbstractFloat}
    # Urban landscape interception
    Rain = Rain + Runon_tm1  # [mm/dth]
    In_imp = In_imp_tm1 + Rain - E_imp * dth * 3600 * 1000 / row  # [mm]
    Lk_imp = min(K_imp, In_imp/dth)  # [mm/h] Leakage from impervious area
    In_imp = In_imp_tm1 + Rain - E_imp * dth * 3600 * 1000 / row - Lk_imp * dth  # [mm]
    q_runon_imp = (In_imp > In_max_imp) * (In_imp - In_max_imp)  # [mm/dth]
    In_imp = In_imp - (In_imp > In_max_imp) * (In_imp - In_max_imp)  # [mm]
    dIn_imp_dt = In_imp - In_imp_tm1  # [mm/dth]

    # Volume Balance check
    WBalance_In_imp =
        Rain - E_imp * dth * 3600 * 1000 / row - q_runon_imp - Lk_imp * dth - dIn_imp_dt  # [mm/dth]

    return q_runon_imp, In_imp, dIn_imp_dt, Lk_imp, WBalance_In_imp
end
