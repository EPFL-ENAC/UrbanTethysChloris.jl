abstract type AbstractBuildingEnergyModelVariables{FT<:AbstractFloat} <:
              AbstractModelVariables{FT} end

"""
    BuildingEnergyModelVariables{FT<:AbstractFloat, N} <: AbstractModelVariables{FT}

Variables for the building energy model.

# Fields
## Building Interior Temperatures (TempVecB)
- `Tceiling`: Building interior ceiling temperature [K]
- `Tinwallsun`: Building interior sunlit wall temperature [K]
- `Tinwallshd`: Building interior shaded wall temperature [K]
- `Twindows`: Building window temperature [K]
- `Tinground`: Building interior ground/floor temperature [K]
- `Tintmass`: Building interior internal heat storage element temperature [K]
- `Tbin`: Building interior air temperature [K]
- `qbin`: Building interior specific humidity [kg/kg]

## Building Interior Humidity (HumidityBuilding)
- `esatbin`: Saturation vapor pressure at building interior temperature [Pa]
- `ebin`: Vapor pressure in building interior [Pa]
- `RHbin`: Relative humidity in building interior [-]

## Building Interior Sensible Heat (HbuildInt)
- `HBinRoof`: Sensible heat flux building interior roof area [W/m² roof area]
- `HbinWallSun`: Sensible heat flux building interior sunlit wall area [W/m² wall area]
- `HbinWallshd`: Sensible heat flux building interior shaded wall area [W/m² wall area]
- `HBinGround`: Sensible heat flux building interior ground area [W/m² ground area]
- `HbinIntMass`: Sensible heat flux building interior mass area [W/m² wall area]
- `HbuildInSurf`: Sensible heat flux building interior total (walls+roof+ceiling+mass) [W/m² ground area]
- `Hvent`: Sensible heat flux building due to ventilation [W/m² ground area]
- `Hequip`: Sensible heat flux building due to equipment [W/m² ground area]
- `Hpeople`: Sensible heat flux building due to people [W/m² ground area]
- `H_AC_Heat`: Sensible heat flux building due to HVAC [W/m² ground area]
- `dSH_air`: Sensible heat flux building due to change in heat storage in air [W/m² ground area]

## Building Interior Latent Heat (LEbuildInt)
- `LEvent`: Latent heat flux building due to ventilation [W/m² ground area]
- `LEequip`: Latent heat flux building due to equipment [W/m² ground area]
- `LEpeople`: Latent heat flux building due to people [W/m² ground area]
- `LE_AC_Heat`: Latent heat flux building due to HVAC [W/m² ground area]
- `dSLE_air`: Latent heat flux building due to change in moisture in air [W/m² ground area]

## Building Interior Conductive Heat Fluxes (GbuildInt)
- `G2Roof`: Conductive heat flux reaching building roof interior [W/m² roof area]
- `G2WallSun`: Conductive heat flux reaching building sunlit wall interior [W/m² wall area]
- `G2WallShade`: Conductive heat flux reaching building shaded wall interior [W/m² wall area]
- `Gfloor`: Conductive heat flux from building floor [W/m² ground area]
- `dSinternalMass`: Change in heat storage in internal mass [W/m² wall area]

## Absorbed Shortwave Radiation in Building (SWRabsB)
- `SWRabsCeiling`: Absorbed shortwave radiation by building interior ceiling [W/m² roof area]
- `SWRabsWallsun`: Absorbed shortwave radiation by building interior sunlit wall [W/m² wall area]
- `SWRabsWallshd`: Absorbed shortwave radiation by building interior shaded wall [W/m² wall area]
- `SWRabsGround`: Absorbed shortwave radiation by building interior ground [W/m² ground area]
- `SWRabsInternalMass`: Absorbed shortwave radiation by building internal mass [W/m² wall area]

## Absorbed Longwave Radiation in Building (LWRabsB)
- `LWRabsCeiling`: Absorbed longwave radiation by building interior ceiling [W/m² roof area]
- `LWRabsWallsun`: Absorbed longwave radiation by building interior sunlit wall [W/m² wall area]
- `LWRabsWallshd`: Absorbed longwave radiation by building interior shaded wall [W/m² wall area]
- `LWRabsGround`: Absorbed longwave radiation by building interior ground [W/m² ground area]
- `LWRabsInternalMass`: Absorbed longwave radiation by building internal mass [W/m² wall area]

## Building Energy Model Waste Heat (BEMWasteHeat)
- `SensibleFromAC_Can`: Sensible heat added to canyon air due to air conditioning energy use [W/m² canyon ground]
- `LatentFromAC_Can`: Latent heat added to canyon air due to air conditioning energy use [W/m² canyon ground]
- `WaterFromAC_Can`: Water that is condensed and removed as runoff in sewer system [W/m² canyon ground]
- `SensibleFromHeat_Can`: Sensible heat added to canyon air due to heating [W/m² canyon ground]
- `LatentFromHeat_Can`: Latent heat added to canyon air due to heating [W/m² canyon ground]
- `SensibleFromVent_Can`: Sensible heat removed or added to the canyon due to exchange of indoor to outdoor air [W/m² canyon ground]
- `LatentFromVent_Can`: Latent heat removed or added to the canyon due to exchange of indoor to outdoor air [W/m² canyon ground]
- `TotAnthInput_URB`: Total anthropogenic heat output to the urban area due to HVAC [W/m² urban]

## Building Energy Use (BEMEnergyUse)
- `EnergyForAC`: Energy consumption for AC [total building interior]
- `EnergyForAC_H`: Energy consumption due to sensible heat load for AC [total building interior]
- `EnergyForAC_LE`: Energy consumption for AC latent [total building interior]
- `EnergyForHeating`: Energy consumption for heating [total building interior]

## AC Parameters Time Series (ParACHeat_ts)
- `AC_on`: Indicating the timesteps in which AC is switched on
- `AC_onCool`: Indicating the timesteps in which AC is switched on due to cooling
- `AC_onDehum`: Indicating the timesteps in which AC is switched on due to dehumidification
- `Heat_on`: Indicating the timesteps in which heating is switched on
"""
Base.@kwdef struct BuildingEnergyModelVariables{FT<:AbstractFloat,N} <:
                   AbstractBuildingEnergyModelVariables{FT}
    # Building Interior Temperatures (TempVecB)
    Tceiling::Array{FT,N}
    Tinwallsun::Array{FT,N}
    Tinwallshd::Array{FT,N}
    Twindows::Array{FT,N}
    Tinground::Array{FT,N}
    Tintmass::Array{FT,N}
    Tbin::Array{FT,N}
    qbin::Array{FT,N}

    # Building Interior Humidity (HumidityBuilding)
    esatbin::Array{FT,N}
    ebin::Array{FT,N}
    RHbin::Array{FT,N}

    # Building Interior Sensible Heat (HbuildInt)
    HBinRoof::Array{FT,N}
    HbinWallSun::Array{FT,N}
    HbinWallshd::Array{FT,N}
    HBinGround::Array{FT,N}
    HbinIntMass::Array{FT,N}
    HbuildInSurf::Array{FT,N}
    Hvent::Array{FT,N}
    Hequip::Array{FT,N}
    Hpeople::Array{FT,N}
    H_AC_Heat::Array{FT,N}
    dSH_air::Array{FT,N}

    # Building Interior Latent Heat (LEbuildInt)
    LEvent::Array{FT,N}
    LEequip::Array{FT,N}
    LEpeople::Array{FT,N}
    LE_AC_Heat::Array{FT,N}
    dSLE_air::Array{FT,N}

    # Building Interior Conductive Heat Fluxes (GbuildInt)
    G2Roof::Array{FT,N}
    G2WallSun::Array{FT,N}
    G2WallShade::Array{FT,N}
    Gfloor::Array{FT,N}
    dSinternalMass::Array{FT,N}

    # Absorbed Shortwave Radiation in Building (SWRabsB)
    SWRabsCeiling::Array{FT,N}
    SWRabsWallsun::Array{FT,N}
    SWRabsWallshd::Array{FT,N}
    SWRabsGround::Array{FT,N}
    SWRabsInternalMass::Array{FT,N}

    # Absorbed Longwave Radiation in Building (LWRabsB)
    LWRabsCeiling::Array{FT,N}
    LWRabsWallsun::Array{FT,N}
    LWRabsWallshd::Array{FT,N}
    LWRabsGround::Array{FT,N}
    LWRabsInternalMass::Array{FT,N}

    # Building Energy Model Waste Heat (BEMWasteHeat)
    SensibleFromAC_Can::Array{FT,N}
    LatentFromAC_Can::Array{FT,N}
    WaterFromAC_Can::Array{FT,N}
    SensibleFromHeat_Can::Array{FT,N}
    LatentFromHeat_Can::Array{FT,N}
    SensibleFromVent_Can::Array{FT,N}
    LatentFromVent_Can::Array{FT,N}
    TotAnthInput_URB::Array{FT,N}

    # Building Energy Use (BEMEnergyUse)
    EnergyForAC::Array{FT,N}
    EnergyForAC_H::Array{FT,N}
    EnergyForAC_LE::Array{FT,N}
    EnergyForHeating::Array{FT,N}

    # AC Parameters Time Series (ParACHeat_ts)
    AC_on::Array{FT,N}
    AC_onCool::Array{FT,N}
    AC_onDehum::Array{FT,N}
    Heat_on::Array{FT,N}
end

function Base.getproperty(
    obj::BuildingEnergyModelVariables{FT,0}, field::Symbol
) where {FT<:AbstractFloat}
    return getfield(obj, field)[]
end

function initialize_building_energy_model_variables(
    ::Type{FT}, N::Int, hours::Int=1, Tatm::FT=FT(293.15), AtmSpecific::FT=FT(0.0)
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        BuildingEnergyModelVariables,
        Dict{String,Any}(),
        (FT, N),
        hours,
        Tatm,
        AtmSpecific,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    Tatm::FT,
    AtmSpecific::FT,
) where {FT<:AbstractFloat,T<:AbstractBuildingEnergyModelVariables}
    processed = Dict{String,Any}()

    dimensions = get_dimensions(T, data, params, hours)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    # Temperature and humidity fields that should be initialized with non-zero values
    temp_fields = [
        "Tceiling", "Tinwallsun", "Tinwallshd", "Twindows", "Tinground", "Tintmass", "Tbin"
    ]
    humidity_fields = ["qbin"]

    if params[2] == 1
        for var in temp_fields
            processed[var][1] = Tatm
        end
        for var in humidity_fields
            processed[var][1] = AtmSpecific
        end
    end

    return processed
end
