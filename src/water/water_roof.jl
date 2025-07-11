"""
    water_roof(
        Eroof_imp::FT,
        Eroof_veg::FT,
        Eroof_ground::FT,
        Eroof_soil::FT,
        TEroof_veg::FT,
        MeteoData::NamedTuple,
        Int_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        Runon_ittm::NamedTuple,
        FractionsRoof::NamedTuple,
        ParSoilRoof::NamedTuple,
        ParCalculation::NamedTuple,
        ParVegRoof::NamedTuple,
        Anthropogenic::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate water balance for roof surfaces including vegetation, impervious areas, and soil.

# Arguments
- `Eroof_imp`: Evaporation from impervious surfaces [kg/m²s]
- `Eroof_veg`: Evaporation from vegetation [kg/m²s]
- `Eroof_ground`: Evaporation from ground under vegetation [kg/m²s]
- `Eroof_soil`: Evaporation from soil [kg/m²s]
- `TEroof_veg`: Transpiration from vegetation [kg/m²s]
- `MeteoData`: Meteorological data
- `Int_ittm`: Previous timestep interception values
- `Owater_ittm`: Previous timestep soil moisture values
- `Runon_ittm`: Previous timestep runon values
- `FractionsRoof`: Roof surface fractions
- `ParSoilRoof`: Soil parameters for roof
- `ParCalculation`: Calculation parameters
- `ParVegRoof`: Vegetation parameters for roof
- `Anthropogenic`: Anthropogenic water inputs

# Returns
- `q_runon_imp`: Runoff from impervious surfaces [mm/dth]
- `In_imp`: Interception on impervious surfaces [mm]
- `dIn_imp_dt`: Change in impervious surface interception [mm/dth]
- `Lk_imp`: Leakage from impervious surfaces [mm/h]
- `q_runon_veg`: Runoff from vegetation [mm/dth]
- `In_veg`: Vegetation interception [mm]
- `dIn_veg_dt`: Change in vegetation interception [mm/dth]
- `q_runon_ground`: Ground runoff [mm/dth]
- `In_ground`: Ground interception [mm]
- `dIn_ground_dt`: Change in ground interception [mm/dth]
- `dV_dt`: Change in soil water volume [mm/dth]
- `f_ground`: Ground infiltration rate [mm/h]
- `V`: Soil water volume [mm]
- `O`: Soil moisture [-]
- `OS`: Surface soil moisture [-]
- `Lk`: Soil leakage [mm/h]
- `Psi_s`: Soil water potential [MPa]
- `Exwat`: Extractable water [mm/h]
- `Rd`: Surface runoff [mm/dth]
- `TEroof_veg`: Updated vegetation transpiration [kg/m²s]
- `Eroof_soil`: Updated soil evaporation [kg/m²s]
- `Runoff`: Total roof runoff [mm/dth]
- `Runon`: Total roof runon [mm/dth]
- `WBalance_In_imp`: Water balance for impervious surfaces [mm/dth]
- `WBalance_In_veg`: Water balance for vegetation [mm/dth]
- `WBalance_In_ground`: Water balance for ground [mm/dth]
- `WBalance_soil`: Water balance for soil [mm/dth]
- `WBalance_imp_tot`: Total water balance for impervious surfaces [mm/dth]
- `WBalance_veg_tot`: Total water balance for vegetation [mm/dth]
- `WBalance_tot`: Total water balance [mm/dth]
"""
function water_roof(
    Eroof_imp::FT,
    Eroof_veg::FT,
    Eroof_ground::FT,
    Eroof_soil::FT,
    TEroof_veg::FT,
    MeteoData::NamedTuple,
    Int_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    Runon_ittm::NamedTuple,
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions,
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParCalculation::NamedTuple,
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    Anthropogenic::NamedTuple,
) where {FT<:AbstractFloat}
    # Extract parameters from NamedTuples
    Rain = MeteoData.Rain
    In_imp_tm1 = Int_ittm.IntRoofImp
    In_veg_tm1 = Int_ittm.IntRoofVegPlant
    In_ground_tm1 = Int_ittm.IntRoofVegGround
    Otm1 = Owater_ittm.OwRoofSoilVeg
    Runon_tm1 = Runon_ittm.RunonRoofTot

    Per_runoff = FractionsRoof.Per_runoff
    fveg = FractionsRoof.fveg
    fimp = FractionsRoof.fimp

    In_max_imp = ParSoilRoof.In_max_imp
    In_max_ground = ParSoilRoof.In_max_ground
    K_imp = ParSoilRoof.Kimp
    Sp_In = ParSoilRoof.Sp_In
    Zs = ParSoilRoof.Zs
    Pcla = ParSoilRoof.Pcla
    Psan = ParSoilRoof.Psan
    Porg = ParSoilRoof.Porg
    Kfc = ParSoilRoof.Kfc
    Phy = ParSoilRoof.Phy
    SPAR = ParSoilRoof.SPAR
    Kbot = ParSoilRoof.Kbot

    dth = ParCalculation.dth
    row = ParCalculation.row

    LAI = ParVegRoof.LAI
    SAI = ParVegRoof.SAI
    Rrootl = ParVegRoof.Rrootl
    PsiL50 = ParVegRoof.PsiL50
    PsiX50 = ParVegRoof.PsiX50
    CASE_ROOT = ParVegRoof.CASE_ROOT
    ZR95 = ParVegRoof.ZR95
    ZR50 = ParVegRoof.ZR50
    ZRmax = ParVegRoof.ZRmax

    # Calculate water fluxes
    q_runon_imp, In_imp, dIn_imp_dt, Lk_imp, WBalance_In_imp = water_impervious(
        Rain, Runon_tm1, Eroof_imp, In_imp_tm1, dth, row, In_max_imp, K_imp
    )

    q_runon_veg, In_veg, dIn_veg_dt, WBalance_In_veg = water_vegetation(
        Rain, Eroof_veg, In_veg_tm1, Sp_In, LAI, SAI, row, dth
    )

    q_runon_ground, In_ground, dIn_ground_dt, f_ground, WBalance_In_ground = water_ground(
        q_runon_veg + Anthropogenic.Waterf_roof,
        Runon_tm1,
        Eroof_ground,
        Otm1,
        In_ground_tm1,
        In_max_ground,
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT,
        CASE_ROOT,
        [0.0],
        ZR95,
        [0.0],
        ZR50,
        [0.0],
        ZRmax,
        Zs,
        dth,
        row,
    )

    V, O, OS, Lk, _, Psi_s, _, Exwat, Rd, TEroof_veg, _, Eroof_soil, dV_dt, WBalance_soil, Psi_soil, Ko = water_soil(
        Otm1,
        f_ground,
        zero(FT),
        TEroof_veg,
        Eroof_soil,
        zeros(FT, length(Zs)),
        dth,
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT,
        CASE_ROOT,
        [0.0],
        ZR95,
        [0.0],
        ZR50,
        [0.0],
        ZRmax,
        [0.0],
        Rrootl,
        [0.0],
        PsiL50,
        [0.0],
        PsiX50,
        Zs,
        row,
    )

    # Calculate runoff and runon
    Runoff = Per_runoff * (fimp * q_runon_imp + fveg * (q_runon_ground + Rd))  # [mm/dth]
    Runon = (1 - Per_runoff) * (fimp * q_runon_imp + fveg * (q_runon_ground + Rd))  # [mm/dth]

    # Water balance checks
    WBalance_imp_tot =
        Rain + Runon_tm1 - Eroof_imp * dth * 3600 * 1000 / row - q_runon_imp -
        Lk_imp * dth - dIn_imp_dt

    WBalance_veg_tot =
        Rain + Runon_tm1 + Anthropogenic.Waterf_roof -
        (Eroof_veg + Eroof_ground + Eroof_soil + TEroof_veg) * dth * 3600 * 1000 / row -
        Lk * dth - q_runon_ground - Rd - dIn_veg_dt - dIn_ground_dt - dV_dt

    E_tot =
        (fimp * Eroof_imp + fveg * (Eroof_veg + Eroof_ground + Eroof_soil + TEroof_veg)) *
        3600 *
        1000 / row  # [mm/h]
    Leak_tot = fimp * Lk_imp + fveg * Lk  # [mm/h]
    Storage_tot = fimp * dIn_imp_dt + fveg * (dIn_veg_dt + dIn_ground_dt + dV_dt)  # [mm/dth]

    WBalance_tot =
        Rain + Runon_tm1 + fveg * Anthropogenic.Waterf_roof - E_tot * dth - Leak_tot * dth -
        Runoff - Runon - Storage_tot

    return q_runon_imp,
    In_imp,
    dIn_imp_dt,
    Lk_imp,
    q_runon_veg,
    In_veg,
    dIn_veg_dt,
    q_runon_ground,
    In_ground,
    dIn_ground_dt,
    dV_dt,
    f_ground,
    V,
    O,
    OS,
    Lk,
    Psi_s,
    Exwat,
    Rd,
    TEroof_veg,
    Eroof_soil,
    Runoff,
    Runon,
    WBalance_In_imp,
    WBalance_In_veg,
    WBalance_In_ground,
    WBalance_soil,
    WBalance_imp_tot,
    WBalance_veg_tot,
    WBalance_tot
end
