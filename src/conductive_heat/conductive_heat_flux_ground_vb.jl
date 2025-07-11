"""
    conductive_heat_flux_ground_vb(
        TemperatureC::Vector{FT},
        TempDamp_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        TempVec_ittm::NamedTuple,
        ParCalculation::NamedTuple,
        ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        type::Int
    ) where {FT<:AbstractFloat}

Calculate conductive heat flux for bare/vegetated ground fractions.

# Arguments
- `TemperatureC`: Canyon temperatures vector
- `TempDamp_ittm`: Temperature damping vectors at previous time step
- `Owater_ittm`: Water content vectors at previous time step
- `TempVec_ittm`: Temperature vectors at previous time step
- `ParCalculation`: Calculation parameters
- `ParSoilGround`: Ground soil parameters
- `ParVegGround`: Ground vegetation parameters
- `ParVegTree`: Tree vegetation parameters
- `FractionsGround`: Ground surface fractions
- `type`: 0 for bare ground, 1 for vegetated ground

# Returns
- `G::FT`: Ground heat flux [W/m²]
- `Tdp::FT`: Damping temperature [K]
"""
function conductive_heat_flux_ground_vb(
    TemperatureC::Vector{FT},
    TempDamp_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    TempVec_ittm::NamedTuple,
    ParCalculation::NamedTuple,
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    type::Int,
) where {FT<:AbstractFloat}
    if type == 0  # bare soil
        Ts = TemperatureC[2]
        Tdptm1 = TempDamp_ittm.TDampGroundBare
        Otm1 = Owater_ittm.OwGroundSoilBare
        Tstm1 = TempVec_ittm.TGroundBare
        Csoil = FractionsGround.fbare > 0
    elseif type == 1  # vegetated ground
        Ts = TemperatureC[3]
        Tdptm1 = TempDamp_ittm.TDampGroundVeg
        Otm1 = Owater_ittm.OwGroundSoilVeg
        Tstm1 = TempVec_ittm.TGroundVeg
        Csoil = FractionsGround.fveg > 0
    else
        throw(
            ArgumentError(
                "Please enter a valid specification of ground type. 0 = bare, 1 = vegetated"
            ),
        )
    end

    dts = ParCalculation.dts
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
        ZR95_H,
        ZR95_L,
        ZR50_H,
        ZR50_L,
        ZRmax_H,
        ZRmax_L,
        Zs,
    )

    Otm1Ave = sum((dz ./ sum(dz)) .* Otm1)

    _, _, CTt = soil_thermal_properties(
        Tdptm1 - 273.15, rsd, lan_dry, lan_s, cv_s, Osat, Ohy, [Otm1Ave]
    )

    G, Tdp = soil_heat(dts, Ts - 273.15, Tstm1 - 273.15, Tdptm1 - 273.15, CTt[1])
    Tdp = Tdp + 273.15
    G = G * Csoil

    return G, Tdp
end
