"""
    IndoorOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}

Optical properties for indoor building surfaces.

# Fields
- `abc::FT`: Albedo ceiling (-)
- `abw::FT`: Albedo wall (-)
- `abg::FT`: Albedo ground (-)
- `abm::FT`: Albedo internal mass (-)
- `ec::FT`: Emissivity ceiling (-)
- `eg::FT`: Emissivity ground (-)
- `ew::FT`: Emissivity wall (-)
- `em::FT`: Emissivity internal mass (-)
"""
Base.@kwdef struct IndoorOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}
    abc::FT
    abw::FT
    abg::FT
    abm::FT
    ec::FT
    eg::FT
    ew::FT
    em::FT
end

function initialize_indooropticalproperties(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, IndoorOpticalProperties, data)
end

function IndoorOpticalProperties(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return IndoorOpticalProperties{FT}(
        data["abc"],
        data["abw"],
        data["abg"],
        data["abm"],
        data["ec"],
        data["eg"],
        data["ew"],
        data["em"],
    )
end

"""
    ThermalBuilding{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters specifying thermal properties and dimensions of building internal mass.

# Fields
- `IntMassOn::Int`: Include building internal mass in calculation (0/1)
- `FloorHeight::FT`: Average floor height in building (m)
- `dzFloor::FT`: Average thickness of floors in building (m)
- `dzWall::FT`: Average thickness of walls in building (m)
- `lan_ground_floor::FT`: Building ground thermal conductivity (W/m K)
- `cv_ground_floor::FT`: Building ground volumetric heat capacity (J/m³ K)
- `lan_floor_IntMass::FT`: Internal mass floor thermal conductivity (W/m K)
- `cv_floor_IntMass::FT`: Internal mass floor volumetric heat capacity (J/m³ K)
- `lan_wall_IntMass::FT`: Internal mass wall thermal conductivity (W/m K)
- `cv_wall_IntMass::FT`: Internal mass wall volumetric heat capacity (J/m³ K)
"""
Base.@kwdef struct ThermalBuilding{FT<:AbstractFloat} <: AbstractParameters{FT}
    IntMassOn::Bool
    FloorHeight::FT
    dzFloor::FT
    dzWall::FT
    lan_ground_floor::FT
    cv_ground_floor::FT
    lan_floor_IntMass::FT
    cv_floor_IntMass::FT
    lan_wall_IntMass::FT
    cv_wall_IntMass::FT
end

function initialize_thermalbuildingparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, ThermalBuilding, data)
end

function ThermalBuilding(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return ThermalBuilding{FT}(
        data["IntMassOn"],
        data["FloorHeight"],
        data["dzFloor"],
        data["dzWall"],
        data["lan_ground_floor"],
        data["cv_ground_floor"],
        data["lan_floor_IntMass"],
        data["cv_floor_IntMass"],
        data["lan_wall_IntMass"],
        data["cv_wall_IntMass"],
    )
end

"""
    WindowParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for building windows.

# Fields
- `WindowsOn::Int`: Include windows in simulation (0/1)
- `GlazingRatio::FT`: Window-to-wall ratio (-)
- `Uvalue::FT`: U-value of windows (W/m² K)
- `lan_windows::FT`: Thermal conductivity of windows (W/m K)
- `cv_glass::FT`: Volumetric heat capacity of glass (J/m³ K)
- `dztot::FT`: Total thickness of glass layers (m)
- `SHGC::FT`: Solar heat gain coefficient (-)
- `SolarTransmittance::FT`: Solar radiation transmittance through windows (-)
- `SolarAbsorptivity::FT`: Fraction of solar radiation absorbed by window (-)
- `SolarAlbedo::FT`: Window albedo (-)
"""
Base.@kwdef struct WindowParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    WindowsOn::Int
    GlazingRatio::FT
    Uvalue::FT
    lan_windows::FT
    cv_glass::FT
    dztot::FT
    SHGC::FT
    SolarTransmittance::FT
    SolarAbsorptivity::FT
    SolarAlbedo::FT
end

function initialize_windowparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, WindowParameters, data)
end

function WindowParameters(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return WindowParameters{FT}(
        Int(data["WindowsOn"]),
        data["GlazingRatio"],
        data["Uvalue"],
        data["lan_windows"],
        data["cv_glass"],
        data["dztot"],
        data["SHGC"],
        data["SolarTransmittance"],
        data["SolarAbsorptivity"],
        data["SolarAlbedo"],
    )
end

"""
    HVACParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for HVAC system.

# Fields
- `ACon::Bool`: Enable air conditioning (true/false)
- `AC_onCool::Bool`: Enable cooling mode (true/false)
- `AC_onDehum::Bool`: Enable dehumidification mode (true/false)
- `MasterOn::Bool`: Enable master control (true/false)
- `Heatingon::Bool`: Enable heating (true/false)
- `TsetpointCooling::FT`: Cooling setpoint temperature (K)
- `TsetpointHeating::FT`: Heating setpoint temperature (K)
- `RHsetpointCooling::FT`: Cooling setpoint relative humidity (%)
- `RHsetpointHeating::FT`: Heating setpoint relative humidity (%)
- `ACH::FT`: Air changes per hour (1/h)
- `COPAC::FT`: Coefficient of performance for AC (-)
- `COPHeat::FT`: Coefficient of performance for heating (-)
- `f_ACLatentToQ::FT`: Fraction of latent heat removed by AC that is condensed (-)
- `q_RHspCooling::FT`: Latent heat for cooling relative humidity setpoint (W/m²)
"""
Base.@kwdef struct HVACParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    ACon::Bool
    AC_onCool::Bool=false
    AC_onDehum::Bool=false
    MasterOn::Bool=false
    Heatingon::Bool
    TsetpointCooling::FT
    TsetpointHeating::FT
    RHsetpointCooling::FT
    RHsetpointHeating::FT
    ACH::FT
    COPAC::FT
    COPHeat::FT
    f_ACLatentToQ::FT
    q_RHspCooling::FT=zero(FT)
end

function TethysChlorisCore.get_optional_fields(::Type{HVACParameters})
    return [:AC_onCool, :AC_onDehum, :MasterOn, :q_RHspCooling]
end

function initialize_hvacparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, HVACParameters, data)
end

function HVACParameters(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return HVACParameters{FT}(
        data["ACon"],
        data["AC_onCool"],
        data["AC_onDehum"],
        data["MasterOn"],
        data["Heatingon"],
        data["TsetpointCooling"],
        data["TsetpointHeating"],
        data["RHsetpointCooling"],
        data["RHsetpointHeating"],
        data["ACH"],
        data["COPAC"],
        data["COPHeat"],
        data["f_ACLatentToQ"],
        data["q_RHspCooling"],
    )
end

"""
    BuildingEnergyModelParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Container for all building energy model parameters.

# Fields
- `indoor_optical::IndoorOpticalProperties{FT}`: Indoor surface optical properties
- `thermal::ThermalBuilding{FT}`: Building thermal properties
- `windows::WindowParameters{FT}`: Window parameters
- `hvac::HVACParameters{FT}`: HVAC system parameters
"""
Base.@kwdef struct BuildingEnergyModelParameters{FT<:AbstractFloat} <:
                   AbstractParameters{FT}
    indoor_optical::IndoorOpticalProperties{FT}
    thermal::ThermalBuilding{FT}
    windows::WindowParameters{FT}
    hvac::HVACParameters{FT}
end

function initialize_building_energy_model_parameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["indoor_optical"] = initialize_indooropticalproperties(
        FT, data["indoor_optical"]
    )
    processed["thermal"] = initialize_thermalbuildingparameters(FT, data["thermal"])
    processed["windows"] = initialize_windowparameters(FT, data["windows"])
    processed["hvac"] = initialize_hvacparameters(FT, data["hvac"])

    return initialize(FT, BuildingEnergyModelParameters, processed, (FT,))
end
