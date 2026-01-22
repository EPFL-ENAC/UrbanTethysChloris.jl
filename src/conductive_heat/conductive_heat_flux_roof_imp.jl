"""
    conductive_heat_flux_roof_imp(
        TemperatureR::Vector{FT},
        TemperatureB::Union{FT, Vector{FT}},
        TempVec_ittm::NamedTuple,
        Anthropogenic::NamedTuple,
        ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParCalculation::NamedTuple,
        BEM_on::Bool
    ) where {FT<:AbstractFloat}

Calculate conductive heat flux through an impervious roof.

# Arguments
- `TemperatureR`: Surface temperature of roof [K]
- `TemperatureB`: Interior temperature in the building [K]
- `TempVec_ittm`: Temperature vectors at previous time step
- `Anthropogenic`: Anthropogenic parameters
- `ParThermalRoof`: Thermal parameters for roof
- `ParSoilRoof`: Soil parameters for roof
- `ParCalculation`: Calculation parameters
- `BEM_on`: Building Energy Model switch (true/false)

# Returns
- `G1`: Heat flux from surface to concrete interior [W/m²]
- `G2`: Heat flux from concrete interior to building interior [W/m²]
- `dS`: Energy storage in the roof [J/m²]
"""
function conductive_heat_flux_roof_imp(
    TemperatureR::Vector{FT},
    TemperatureB::Union{FT,Vector{FT}},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    Anthropogenic::ModelComponents.ForcingInputs.AnthropogenicInputs{FT,0},
    ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParCalculation::NamedTuple,
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Extract inputs
    Ts = TemperatureR[1]
    Tint = TemperatureR[3]
    Tb = BEM_on ? TemperatureB[1] : Anthropogenic.Tb
    Tint_tm1 = TempVec_ittm.TRoofIntImp
    lan_dry1 = ParThermalRoof.lan_dry
    lan_dry2 = ParThermalRoof.lan_dry
    dz1 = ParSoilRoof.dz1
    dz2 = ParSoilRoof.dz2
    cv_s1 = ParThermalRoof.cv_s
    cv_s2 = ParThermalRoof.cv_s
    dts = ParCalculation.dts

    # Computation of heat fluxes
    G1 = lan_dry1 * (Ts - Tint) / dz1
    G2 = lan_dry2 * (Tint - Tb) / dz2
    dS = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

    return G1, G2, dS
end
function conductive_heat_flux_roof_imp(
    TemperatureR::Vector{FT},
    TemperatureB::Union{FT,Vector{FT}},
    TempVec_ittm::NamedTuple,
    Anthropogenic::NamedTuple,
    ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParCalculation::NamedTuple,
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Extract inputs
    Ts = TemperatureR[1]
    Tint = TemperatureR[3]
    Tb = BEM_on ? TemperatureB[1] : Anthropogenic.Tb
    Tint_tm1 = TempVec_ittm.TRoofIntImp
    lan_dry1 = ParThermalRoof.lan_dry
    lan_dry2 = ParThermalRoof.lan_dry
    dz1 = ParSoilRoof.dz1
    dz2 = ParSoilRoof.dz2
    cv_s1 = ParThermalRoof.cv_s
    cv_s2 = ParThermalRoof.cv_s
    dts = ParCalculation.dts

    # Computation of heat fluxes
    G1 = lan_dry1 * (Ts - Tint) / dz1
    G2 = lan_dry2 * (Tint - Tb) / dz2
    dS = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

    return G1, G2, dS
end
