# Subtype for the set of subsets of building energy model variables
abstract type AbstractBuildingEnergyModelVariables{FT<:AbstractFloat} <:
              AbstractModelVariableSet{FT} end

# Subtype for each subset of building energy model variables
abstract type AbstractBuildingEnergyModelVariablesSubset{FT<:AbstractFloat,N} <:
              AbstractModelVariables{FT} end

abstract type AbstractTempVecB{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractHumidityBuilding{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractHbuildInt{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractLEbuildInt{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractGbuildInt{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractSWRabsB{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractLWRabsB{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractBEMWasteHeat{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractBEMEnergyUse{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end
abstract type AbstractParACHeat_ts{FT<:AbstractFloat,N} <:
              AbstractBuildingEnergyModelVariablesSubset{FT,N} end

"""
    TempVecB{FT<:AbstractFloat, N} <: AbstractTempVecB{FT}

Building interior temperatures.

# Fields
- `Tceiling`: Building interior ceiling temperature [K]
- `Tinwallsun`: Building interior sunlit wall temperature [K]
- `Tinwallshd`: Building interior shaded wall temperature [K]
- `Twindows`: Building window temperature [K]
- `Tinground`: Building interior ground/floor temperature [K]
- `Tintmass`: Building interior internal heat storage element temperature [K]
- `Tbin`: Building interior air temperature [K]
- `qbin`: Building interior specific humidity [kg/kg]
"""
Base.@kwdef struct TempVecB{FT<:AbstractFloat,N} <: AbstractTempVecB{FT,N}
    Tceiling::Array{FT,N}
    Tinwallsun::Array{FT,N}
    Tinwallshd::Array{FT,N}
    Twindows::Array{FT,N}
    Tinground::Array{FT,N}
    Tintmass::Array{FT,N}
    Tbin::Array{FT,N}
    qbin::Array{FT,N}
end

"""
    HumidityBuilding{FT<:AbstractFloat, N} <: AbstractHumidityBuilding{FT,N}

Building interior humidity variables.

# Fields
- `esatbin`: Saturation vapor pressure at building interior temperature [Pa]
- `ebin`: Vapor pressure in building interior [Pa]
- `RHbin`: Relative humidity in building interior [-]
"""
Base.@kwdef struct HumidityBuilding{FT<:AbstractFloat,N} <: AbstractHumidityBuilding{FT,N}
    esatbin::Array{FT,N}
    ebin::Array{FT,N}
    RHbin::Array{FT,N}
end

"""
    HbuildInt{FT<:AbstractFloat, N} <: AbstractHbuildInt{FT,N}

Building interior sensible heat fluxes.

# Fields
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
"""
Base.@kwdef struct HbuildInt{FT<:AbstractFloat,N} <: AbstractHbuildInt{FT,N}
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
end

"""
    LEbuildInt{FT<:AbstractFloat, N} <: AbstractLEbuildInt{FT,N}

Building interior latent heat fluxes.

# Fields
- `LEvent`: Latent heat flux building due to ventilation [W/m² ground area]
- `LEequip`: Latent heat flux building due to equipment [W/m² ground area]
- `LEpeople`: Latent heat flux building due to people [W/m² ground area]
- `LE_AC_Heat`: Latent heat flux building due to HVAC [W/m² ground area]
- `dSLE_air`: Latent heat flux building due to change in moisture in air [W/m² ground area]
"""
Base.@kwdef struct LEbuildInt{FT<:AbstractFloat,N} <: AbstractLEbuildInt{FT,N}
    LEvent::Array{FT,N}
    LEequip::Array{FT,N}
    LEpeople::Array{FT,N}
    LE_AC_Heat::Array{FT,N}
    dSLE_air::Array{FT,N}
end

"""
    GbuildInt{FT<:AbstractFloat, N} <: AbstractGbuildInt{FT,N}

Building interior conductive heat fluxes.

# Fields
- `G2Roof`: Conductive heat flux reaching building roof interior [W/m² roof area]
- `G2WallSun`: Conductive heat flux reaching building sunlit wall interior [W/m² wall area]
- `G2WallShade`: Conductive heat flux reaching building shaded wall interior [W/m² wall area]
- `Gfloor`: Conductive heat flux from building floor [W/m² ground area]
- `dSinternalMass`: Change in heat storage in internal mass [W/m² wall area]
"""
Base.@kwdef struct GbuildInt{FT<:AbstractFloat,N} <: AbstractGbuildInt{FT,N}
    G2Roof::Array{FT,N}
    G2WallSun::Array{FT,N}
    G2WallShade::Array{FT,N}
    Gfloor::Array{FT,N}
    dSinternalMass::Array{FT,N}
end

"""
    SWRabsB{FT<:AbstractFloat, N} <: AbstractSWRabsB{FT,N}

Absorbed shortwave radiation in building interior.

# Fields
- `SWRabsCeiling`: Absorbed shortwave radiation by building interior ceiling [W/m² roof area]
- `SWRabsWallsun`: Absorbed shortwave radiation by building interior sunlit wall [W/m² wall area]
- `SWRabsWallshd`: Absorbed shortwave radiation by building interior shaded wall [W/m² wall area]
- `SWRabsGround`: Absorbed shortwave radiation by building interior ground [W/m² ground area]
- `SWRabsInternalMass`: Absorbed shortwave radiation by building internal mass [W/m² wall area]
"""
Base.@kwdef struct SWRabsB{FT<:AbstractFloat,N} <: AbstractSWRabsB{FT,N}
    SWRabsCeiling::Array{FT,N}
    SWRabsWallsun::Array{FT,N}
    SWRabsWallshd::Array{FT,N}
    SWRabsGround::Array{FT,N}
    SWRabsInternalMass::Array{FT,N}
end

"""
    LWRabsB{FT<:AbstractFloat, N} <: AbstractLWRabsB{FT,N}

Absorbed longwave radiation in building interior.

# Fields
- `LWRabsCeiling`: Absorbed longwave radiation by building interior ceiling [W/m² roof area]
- `LWRabsWallsun`: Absorbed longwave radiation by building interior sunlit wall [W/m² wall area]
- `LWRabsWallshd`: Absorbed longwave radiation by building interior shaded wall [W/m² wall area]
- `LWRabsGround`: Absorbed longwave radiation by building interior ground [W/m² ground area]
- `LWRabsInternalMass`: Absorbed longwave radiation by building internal mass [W/m² wall area]
"""
Base.@kwdef struct LWRabsB{FT<:AbstractFloat,N} <: AbstractLWRabsB{FT,N}
    LWRabsCeiling::Array{FT,N}
    LWRabsWallsun::Array{FT,N}
    LWRabsWallshd::Array{FT,N}
    LWRabsGround::Array{FT,N}
    LWRabsInternalMass::Array{FT,N}
end

"""
    BEMWasteHeat{FT<:AbstractFloat, N} <: AbstractBEMWasteHeat{FT,N}

Building energy model waste heat variables.

# Fields
- `SensibleFromAC_Can`: Sensible heat added to canyon air due to air conditioning energy use [W/m² canyon ground]
- `LatentFromAC_Can`: Latent heat added to canyon air due to air conditioning energy use [W/m² canyon ground]
- `WaterFromAC_Can`: Water that is condensed and removed as runoff in sewer system [W/m² canyon ground]
- `SensibleFromHeat_Can`: Sensible heat added to canyon air due to heating [W/m² canyon ground]
- `LatentFromHeat_Can`: Latent heat added to canyon air due to heating [W/m² canyon ground]
- `SensibleFromVent_Can`: Sensible heat removed or added to the canyon due to exchange of indoor to outdoor air [W/m² canyon ground]
- `LatentFromVent_Can`: Latent heat removed or added to the canyon due to exchange of indoor to outdoor air [W/m² canyon ground]
- `TotAnthInput_URB`: Total anthropogenic heat output to the urban area due to HVAC [W/m² urban]
"""
Base.@kwdef struct BEMWasteHeat{FT<:AbstractFloat,N} <: AbstractBEMWasteHeat{FT,N}
    SensibleFromAC_Can::Array{FT,N}
    LatentFromAC_Can::Array{FT,N}
    WaterFromAC_Can::Array{FT,N}
    SensibleFromHeat_Can::Array{FT,N}
    LatentFromHeat_Can::Array{FT,N}
    SensibleFromVent_Can::Array{FT,N}
    LatentFromVent_Can::Array{FT,N}
    TotAnthInput_URB::Array{FT,N}
end

"""
    BEMEnergyUse{FT<:AbstractFloat, N} <: AbstractBEMEnergyUse{FT,N}

Building energy use variables.

# Fields
- `EnergyForAC`: Energy consumption for AC [total building interior]
- `EnergyForAC_H`: Energy consumption due to sensible heat load for AC [total building interior]
- `EnergyForAC_LE`: Energy consumption for AC latent [total building interior]
- `EnergyForHeating`: Energy consumption for heating [total building interior]
"""
Base.@kwdef struct BEMEnergyUse{FT<:AbstractFloat,N} <: AbstractBEMEnergyUse{FT,N}
    EnergyForAC::Array{FT,N}
    EnergyForAC_H::Array{FT,N}
    EnergyForAC_LE::Array{FT,N}
    EnergyForHeating::Array{FT,N}
end

"""
    ParACHeat_ts{FT<:AbstractFloat, N} <: AbstractParACHeat_ts{FT,N}

AC parameters time series variables.

# Fields
- `AC_on`: Indicating the timesteps in which AC is switched on
- `AC_onCool`: Indicating the timesteps in which AC is switched on due to cooling
- `AC_onDehum`: Indicating the timesteps in which AC is switched on due to dehumidification
- `Heat_on`: Indicating the timesteps in which heating is switched on
"""
Base.@kwdef struct ParACHeat_ts{FT<:AbstractFloat,N} <: AbstractParACHeat_ts{FT,N}
    AC_on::Array{FT,N}
    AC_onCool::Array{FT,N}
    AC_onDehum::Array{FT,N}
    Heat_on::Array{FT,N}
end

"""
    BuildingEnergyModelVariables{FT<:AbstractFloat, N} <: AbstractBuildingEnergyModelVariables{FT,N}

Container for all building energy model variable components.

# Fields
- `TempVecB`: Building interior temperatures
- `HumidityBuilding`: Building interior humidity variables
- `HbuildInt`: Building interior sensible heat fluxes
- `LEbuildInt`: Building interior latent heat fluxes
- `GbuildInt`: Building interior conductive heat fluxes
- `SWRabsB`: Absorbed shortwave radiation in building interior
- `LWRabsB`: Absorbed longwave radiation in building interior
- `BEMWasteHeat`: Building energy model waste heat variables
- `BEMEnergyUse`: Building energy use variables
- `ParACHeat_ts`: AC parameters time series variables
"""
Base.@kwdef struct BuildingEnergyModelVariables{FT<:AbstractFloat,N} <:
                   AbstractBuildingEnergyModelVariables{FT}
    TempVecB::TempVecB{FT,N}
    HumidityBuilding::HumidityBuilding{FT,N}
    HbuildInt::HbuildInt{FT,N}
    LEbuildInt::LEbuildInt{FT,N}
    GbuildInt::GbuildInt{FT,N}
    SWRabsB::SWRabsB{FT,N}
    LWRabsB::LWRabsB{FT,N}
    BEMWasteHeat::BEMWasteHeat{FT,N}
    BEMEnergyUse::BEMEnergyUse{FT,N}
    ParACHeat_ts::ParACHeat_ts{FT,N}
end

# Base getproperty methods for scalar access
function Base.getproperty(
    obj::T, field::Symbol
) where {FT<:AbstractFloat,T<:AbstractBuildingEnergyModelVariablesSubset{FT,0}}
    return getfield(obj, field)[]
end

# Initialization functions for individual components
function initialize_tempvecb(
    ::Type{FT}, N::Int, hours::Int=1, Tatm::FT=FT(293.15), AtmSpecific::FT=FT(0.0)
) where {FT<:AbstractFloat}
    return initialize(FT, TempVecB, Dict{String,Any}(), (FT, N), hours, Tatm, AtmSpecific)
end

function initialize_humidity_building(
    ::Type{FT}, N::Int, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, HumidityBuilding, Dict{String,Any}(), (FT, N), hours)
end

function initialize_hbuildint(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, HbuildInt, Dict{String,Any}(), (FT, N), hours)
end

function initialize_lebuildint(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, LEbuildInt, Dict{String,Any}(), (FT, N), hours)
end

function initialize_gbuildint(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, GbuildInt, Dict{String,Any}(), (FT, N), hours)
end

function initialize_swrabsb(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, SWRabsB, Dict{String,Any}(), (FT, N), hours)
end

function initialize_lwrabsb(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, LWRabsB, Dict{String,Any}(), (FT, N), hours)
end

function initialize_bemwasteheat(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, BEMWasteHeat, Dict{String,Any}(), (FT, N), hours)
end

function initialize_bemenergyuse(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, BEMEnergyUse, Dict{String,Any}(), (FT, N), hours)
end

function initialize_paracheat_ts(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, ParACHeat_ts, Dict{String,Any}(), (FT, N), hours)
end

# Main initialization function for BuildingEnergyModelVariables
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

# Preprocess methods for the individual components
function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{TempVecB},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    Tatm::FT,
    AtmSpecific::FT,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(TempVecB, data, params, hours)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    # Initialize temperature and humidity fields
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

# Main preprocess_fields method for BuildingEnergyModelVariables
function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{BuildingEnergyModelVariables},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    Tatm::FT,
    AtmSpecific::FT,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Initialize each component
    processed["TempVecB"] = initialize_tempvecb(FT, params[2], hours, Tatm, AtmSpecific)
    processed["HumidityBuilding"] = initialize_humidity_building(FT, params[2], hours)
    processed["HbuildInt"] = initialize_hbuildint(FT, params[2], hours)
    processed["LEbuildInt"] = initialize_lebuildint(FT, params[2], hours)
    processed["GbuildInt"] = initialize_gbuildint(FT, params[2], hours)
    processed["SWRabsB"] = initialize_swrabsb(FT, params[2], hours)
    processed["LWRabsB"] = initialize_lwrabsb(FT, params[2], hours)
    processed["BEMWasteHeat"] = initialize_bemwasteheat(FT, params[2], hours)
    processed["BEMEnergyUse"] = initialize_bemenergyuse(FT, params[2], hours)
    processed["ParACHeat_ts"] = initialize_paracheat_ts(FT, params[2], hours)

    return processed
end
