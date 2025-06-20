"""
    conductive_heat_flux_ground_fr(
        TemperatureC::Vector{FT},
        TempDamp_ittm::NamedTuple,
        TempVec_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        ParCalculation::NamedTuple,
        ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT}
    ) where {FT<:AbstractFloat}

Calculate conductive heat flux for ground with fraction-resolved temperature.

# Arguments
- `TemperatureC`: Canyon temperatures vector
- `TempDamp_ittm`: Temperature damping vectors at previous time step
- `TempVec_ittm`: Temperature vectors at previous time step
- `Owater_ittm`: Water content vectors at previous time step
- `ParCalculation`: Calculation parameters
- `ParThermalGround`: Ground thermal parameters
- `FractionsGround`: Ground surface fractions
- `ParSoilGround`: Ground soil parameters
- `ParVegTree`: Tree vegetation parameters
- `ParVegGround`: Ground vegetation parameters

# Returns
- `G::FT`: Ground heat flux [W/m²]
- `Tdp::FT`: Damping temperature [K]
"""
function conductive_heat_flux_ground_fr(
    TemperatureC::Vector{FT},
    TempDamp_ittm::NamedTuple,
    TempVec_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    ParCalculation::NamedTuple,
    ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
) where {FT<:AbstractFloat}
    # Extract inputs
    Ts = TemperatureC[1]
    Tdptm1 = TempDamp_ittm.TDampGroundImp
    Tstm1 = TempVec_ittm.TGroundImp
    dts = ParCalculation.dts
    lan_dry_imp = ParThermalGround.lan_dry
    cv_s_imp = ParThermalGround.cv_s
    Cimp = FractionsGround.fimp > 0

    # To calculate the influence of soil moisture on the soil heat capacity
    Otm1 = Owater_ittm.OwGroundSoilImp
    Pcla = ParSoilGround.Pcla
    Psan = ParSoilGround.Psan
    Porg = ParSoilGround.Porg
    Kfc = ParSoilGround.Kfc
    Phy = ParSoilGround.Phy
    SPAR = ParSoilGround.SPAR
    Kbot = ParSoilGround.Kbot
    CASE_ROOT_H = ParVegTree.CASE_ROOT
    CASE_ROOT_L = ParVegGround.CASE_ROOT
    ZR95_H = ParVegTree.ZR95
    ZR95_L = ParVegGround.ZR95
    ZR50_H = ParVegTree.ZR50
    ZR50_L = ParVegGround.ZR50
    ZRmax_H = ParVegTree.ZRmax
    ZRmax_L = ParVegGround.ZRmax
    Zs = ParSoilGround.Zs

    # Calculate soil properties underneath impervious surface
    _, dz, _, Osat, Ohy, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, rsd, lan_dry, lan_s, cv_s = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT_H,
        CASE_ROOT_L,
        [ZR95_H],
        [ZR95_L],
        [ZR50_H],
        [ZR50_L],
        [ZRmax_H],
        [ZRmax_L],
        Zs,
    )

    NotImp = .!isnan.(Otm1) .&& Otm1 .!= 0
    Otm1Ave = sum((dz[NotImp] ./ sum(dz[NotImp])) .* Otm1[NotImp])

    lanSv, cv_Soilv, _ = soil_thermal_properties(
        Tdptm1 - 273.15, rsd, lan_dry, lan_s, cv_s, Osat, Ohy, [Otm1Ave]
    )
    lanS = lanSv[1]
    cv_Soil = cv_Soilv[1]

    # Calculate total thermal properties
    lan_tot = (lan_dry_imp * sum(dz[.!NotImp]) + lanS * sum(dz[NotImp])) / sum(dz)
    cv_tot = (cv_s_imp * sum(dz[.!NotImp]) + cv_Soil * sum(dz[NotImp])) / sum(dz)

    tau = FT(86400)  # [s] time constant
    CTt = 2 * sqrt(π / (lan_tot * cv_tot * tau))  # [K m²/J] Total Thermal Capacity Soil

    G, Tdp = soil_heat(dts, Ts - 273.15, Tstm1 - 273.15, Tdptm1 - 273.15, CTt)
    Tdp = Tdp + 273.15
    G = G * Cimp

    return G, Tdp
end
