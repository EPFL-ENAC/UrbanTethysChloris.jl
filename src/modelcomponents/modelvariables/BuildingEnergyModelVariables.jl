"""
    TempVecB{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct TempVecB{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    Tceiling::FT
    Tinwallsun::FT
    Tinwallshd::FT
    Twindows::FT
    Tinground::FT
    Tintmass::FT
    Tbin::FT
    qbin::FT
    #TODO: missing Tatm field?
end

function TempVecB(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, TempVecB, Dict{String,Any}())
end

function TempVecB(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return TempVecB{FT}(;
        Tceiling=data["Tceiling"],
        Tinwallsun=data["Tinwallsun"],
        Tinwallshd=data["Tinwallshd"],
        Twindows=data["Twindows"],
        Tinground=data["Tinground"],
        Tintmass=data["Tintmass"],
        Tbin=data["Tbin"],
        qbin=data["qbin"],
    )
end

function update!(x::TempVecB{FT}, Ttot::Vector{FT}) where {FT<:AbstractFloat}
    x.Tceiling = Ttot[15]
    x.Tinwallsun = Ttot[16]
    x.Tinwallshd = Ttot[17]
    x.Twindows = Ttot[18]
    x.Tinground = Ttot[19]
    x.Tintmass = Ttot[20]
    x.Tbin = Ttot[21]
    x.qbin = Ttot[22]
    return nothing
end

"""
    HumidityBuilding{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Building interior humidity variables.

# Fields
- `esatbin`: Saturation vapor pressure at building interior temperature [Pa]
- `ebin`: Vapor pressure in building interior [Pa]
- `RHbin`: Relative humidity in building interior [-]
"""
Base.@kwdef mutable struct HumidityBuilding{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    esatbin::FT
    ebin::FT
    RHbin::FT
end

function HumidityBuilding(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, HumidityBuilding, Dict{String,Any}())
end

"""
    HbuildInt{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct HbuildInt{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    HBinRoof::FT
    HbinWallSun::FT
    HbinWallshd::FT
    HBinGround::FT
    HbinIntMass::FT
    HbuildInSurf::FT
    Hvent::FT
    Hequip::FT
    Hpeople::FT
    H_AC_Heat::FT
    dSH_air::FT
end

function HbuildInt(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, HbuildInt, Dict{String,Any}())
end

"""
    LEbuildInt{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Building interior latent heat fluxes.

# Fields
- `LEvent`: Latent heat flux building due to ventilation [W/m² ground area]
- `LEequip`: Latent heat flux building due to equipment [W/m² ground area]
- `LEpeople`: Latent heat flux building due to people [W/m² ground area]
- `LE_AC_Heat`: Latent heat flux building due to HVAC [W/m² ground area]
- `dSLE_air`: Latent heat flux building due to change in moisture in air [W/m² ground area]
"""
Base.@kwdef mutable struct LEbuildInt{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    LEvent::FT
    LEequip::FT
    LEpeople::FT
    LE_AC_Heat::FT
    dSLE_air::FT
end

function LEbuildInt(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, LEbuildInt, Dict{String,Any}())
end

"""
    GbuildInt{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Building interior conductive heat fluxes.

# Fields
- `G2Roof`: Conductive heat flux reaching building roof interior [W/m² roof area]
- `G2WallSun`: Conductive heat flux reaching building sunlit wall interior [W/m² wall area]
- `G2WallShade`: Conductive heat flux reaching building shaded wall interior [W/m² wall area]
- `Gfloor`: Conductive heat flux from building floor [W/m² ground area]
- `dSinternalMass`: Change in heat storage in internal mass [W/m² wall area]
"""
Base.@kwdef mutable struct GbuildInt{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    G2Roof::FT
    G2WallSun::FT
    G2WallShade::FT
    Gfloor::FT
    dSinternalMass::FT
end

function GbuildInt(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, GbuildInt, Dict{String,Any}())
end

"""
    SWRabsB{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Absorbed shortwave radiation in building interior.

# Fields
- `SWRabsCeiling`: Absorbed shortwave radiation by building interior ceiling [W/m² roof area]
- `SWRabsWallsun`: Absorbed shortwave radiation by building interior sunlit wall [W/m² wall area]
- `SWRabsWallshd`: Absorbed shortwave radiation by building interior shaded wall [W/m² wall area]
- `SWRabsGround`: Absorbed shortwave radiation by building interior ground [W/m² ground area]
- `SWRabsInternalMass`: Absorbed shortwave radiation by building internal mass [W/m² wall area]
"""
Base.@kwdef mutable struct SWRabsB{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    SWRabsCeiling::FT
    SWRabsWallsun::FT
    SWRabsWallshd::FT
    SWRabsGround::FT
    SWRabsInternalMass::FT
end

function SWRabsB(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, SWRabsB, Dict{String,Any}())
end

"""
    LWRabsB{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Absorbed longwave radiation in building interior.

# Fields
- `LWRabsCeiling`: Absorbed longwave radiation by building interior ceiling [W/m² roof area]
- `LWRabsWallsun`: Absorbed longwave radiation by building interior sunlit wall [W/m² wall area]
- `LWRabsWallshd`: Absorbed longwave radiation by building interior shaded wall [W/m² wall area]
- `LWRabsGround`: Absorbed longwave radiation by building interior ground [W/m² ground area]
- `LWRabsInternalMass`: Absorbed longwave radiation by building internal mass [W/m² wall area]
"""
Base.@kwdef mutable struct LWRabsB{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    LWRabsCeiling::FT
    LWRabsWallsun::FT
    LWRabsWallshd::FT
    LWRabsGround::FT
    LWRabsInternalMass::FT
end

function LWRabsB(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, LWRabsB, Dict{String,Any}())
end

"""
    BEMWasteHeat{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct BEMWasteHeat{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    SensibleFromAC_Can::FT
    LatentFromAC_Can::FT
    WaterFromAC_Can::FT
    SensibleFromHeat_Can::FT
    LatentFromHeat_Can::FT
    SensibleFromVent_Can::FT
    LatentFromVent_Can::FT
    TotAnthInput_URB::FT
end

function BEMWasteHeat(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, BEMWasteHeat, Dict{String,Any}())
end

"""
    BEMEnergyUse{FT<:AbstractFloat} <: AbstractBEMEnergyUse{FT}

Building energy use variables.

# Fields
- `EnergyForAC`: Energy consumption for AC [total building interior]
- `EnergyForAC_H`: Energy consumption due to sensible heat load for AC [total building interior]
- `EnergyForAC_LE`: Energy consumption for AC latent [total building interior]
- `EnergyForHeating`: Energy consumption for heating [total building interior]
"""
Base.@kwdef mutable struct BEMEnergyUse{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    EnergyForAC::FT
    EnergyForAC_H::FT
    EnergyForAC_LE::FT
    EnergyForHeating::FT
end

function BEMEnergyUse(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, BEMEnergyUse, Dict{String,Any}())
end

"""
    ParACHeat_ts{FT<:AbstractFloat} <: AbstractParACHeat_ts{FT}

AC parameters time series variables.

# Fields
- `AC_on`: Indicating the timesteps in which AC is switched on
- `AC_onCool`: Indicating the timesteps in which AC is switched on due to cooling
- `AC_onDehum`: Indicating the timesteps in which AC is switched on due to dehumidification
- `Heat_on`: Indicating the timesteps in which heating is switched on
"""
Base.@kwdef mutable struct ParACHeat_ts{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    AC_on::FT
    AC_onCool::FT
    AC_onDehum::FT
    Heat_on::FT
end

function ParACHeat_ts(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, ParACHeat_ts, Dict{String,Any}())
end

"""
    BuildingEnergyModelVariables{FT<:AbstractFloat} <: AbstractBuildingEnergyModelVariables{FT}

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
Base.@kwdef struct BuildingEnergyModelVariables{FT<:AbstractFloat} <:
                   AbstractModelVariableSet{FT}
    TempVecB::TempVecB{FT}
    HumidityBuilding::HumidityBuilding{FT}
    HbuildInt::HbuildInt{FT}
    LEbuildInt::LEbuildInt{FT}
    GbuildInt::GbuildInt{FT}
    SWRabsB::SWRabsB{FT}
    LWRabsB::LWRabsB{FT}
    BEMWasteHeat::BEMWasteHeat{FT}
    BEMEnergyUse::BEMEnergyUse{FT}
    ParACHeat_ts::ParACHeat_ts{FT}
end

function BuildingEnergyModelVariables(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, BuildingEnergyModelVariables, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{BuildingEnergyModelVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Initialize each component
    processed["TempVecB"] = TempVecB(FT)
    processed["HumidityBuilding"] = HumidityBuilding(FT)
    processed["HbuildInt"] = HbuildInt(FT)
    processed["LEbuildInt"] = LEbuildInt(FT)
    processed["GbuildInt"] = GbuildInt(FT)
    processed["SWRabsB"] = SWRabsB(FT)
    processed["LWRabsB"] = LWRabsB(FT)
    processed["BEMWasteHeat"] = BEMWasteHeat(FT)
    processed["BEMEnergyUse"] = BEMEnergyUse(FT)
    processed["ParACHeat_ts"] = ParACHeat_ts(FT)

    return processed
end
