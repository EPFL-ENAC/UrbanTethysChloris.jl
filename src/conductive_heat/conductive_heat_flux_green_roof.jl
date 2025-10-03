"""
    conductive_heat_flux_green_roof(
        TemperatureR::Vector{FT},
        TemperatureB::Vector{FT},
        TempVec_ittm::NamedTuple,
        Anthropogenic::NamedTuple,
        Owater::NamedTuple,
        ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParCalculation::NamedTuple,
        BEM_on::Bool
    ) where {FT<:AbstractFloat}

Calculate conductive heat flux through a green roof.

# Arguments
- `temperature_r`: Surface temperature of roof [K]
- `temperature_b`: Interior temperature in the building [K]
- `temp_vec_ittm`: Temperature vectors at previous time step
- `anthropogenic`: Anthropogenic parameters
- `owater`: Water content parameters
- `par_veg_roof`: Vegetation parameters for roof
- `par_soil_roof`: Soil parameters for roof
- `par_thermal_roof`: Thermal parameters for roof
- `par_calculation`: Calculation parameters
- `bem_on`: Building Energy Model switch (0/1)

# Returns
- `G1::FT`: Heat flux from surface to concrete interior [W/m²]
- `G2::FT`: Heat flux from concrete interior to building interior [W/m²]
- `dS::FT`: Energy storage in the roof [J/m²]
"""
function conductive_heat_flux_green_roof(
    TemperatureR::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::NamedTuple,
    Anthropogenic::NamedTuple,
    Owater::NamedTuple,
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParCalculation::NamedTuple,
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Extract all input parameters
    Troof = TemperatureR[2]
    Tint = TemperatureR[4]
    Tint_tm1 = TempVec_ittm.TRoofIntVeg
    Tb = BEM_on ? TemperatureB[1] : Anthropogenic.Tb
    Otm1 = Owater.OwRoofSoilVeg
    Rrootl = ParVegRoof.Rrootl
    PsiL50 = ParVegRoof.PsiL50
    PsiX50 = ParVegRoof.PsiX50
    CASE_ROOT = ParVegRoof.CASE_ROOT
    ZR95 = ParVegRoof.ZR95
    ZR50 = ParVegRoof.ZR50
    ZRmax = ParVegRoof.ZRmax
    Pcla = ParSoilRoof.Pcla
    Psan = ParSoilRoof.Psan
    Porg = ParSoilRoof.Porg
    Kfc = ParSoilRoof.Kfc
    Phy = ParSoilRoof.Phy
    SPAR = ParSoilRoof.SPAR
    Kbot = ParSoilRoof.Kbot
    Zs = ParSoilRoof.Zs
    dz1 = ParSoilRoof.dz1
    dz2 = ParSoilRoof.dz2
    cv_s2 = ParThermalRoof.cv_s
    lan_dry2 = ParThermalRoof.lan_dry
    dts = ParCalculation.dts

    # Calculate soil parameters
    _, dz, ms, Osat, Ohy, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, rsd, lan_dry, lan_s, cv_s = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT,
        CASE_ROOT,
        [FT(0)],
        ZR95,
        [FT(0)],
        ZR50,
        [FT(0)],
        ZRmax,
        Zs,
    )

    # Temperature calculations
    Tdamptm1 = Tb
    Tdp = fill(Tdamptm1, ms)
    Tdptm1 = copy(Tdp)

    lanS, cv_soil, _ = soil_thermal_properties(
        Tdptm1 .- 273.15, rsd, lan_dry, lan_s, cv_s, Osat, Ohy, Otm1
    )

    # Average soil parameters according to soil layer thickness
    dz = dz ./ sum(dz)
    lanS = dot(dz, lanS)
    cv_soil = dot(dz, cv_soil)

    # Average roof parameters according to roof thickness
    dz_roof = [dz1, dz2]
    cv_roof = [cv_soil, cv_s2]
    dz_roof = dz_roof ./ sum(dz_roof)
    cv_roof = dot(dz_roof, cv_roof)

    # Computation of heat fluxes
    G1 = lanS * (Troof - Tint) / dz1
    G2 = lan_dry2 * (Tint - Tb) / dz2
    dS = cv_roof * (dz1 + dz2) / dts * (Tint - Tint_tm1)

    return G1, G2, dS
end
