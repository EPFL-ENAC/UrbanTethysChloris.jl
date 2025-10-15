"""
    ForcingInputSet{FT<:AbstractFloat} <: AbstractForcingInputSet{FT}

Forcing inputs for the Urban Tethys-Chloris model.

# Fields
- `datetime::Vector{DateTime}`: Date and time of the forcing inputs.
- `anthropogenic::AnthropogenicInputs{FT}`: Anthropogenic inputs.
- `hvacschedule::HVACSchedule{FT}`: HVAC schedule.
- `meteorological::MeteorologicalInputs{FT}`: Meteorological inputs.
- `sunposition::SunPositionInputs{FT}`: Sun position inputs.
"""
Base.@kwdef struct ForcingInputSet{FT<:AbstractFloat} <: AbstractForcingInputSet{FT}
    datetime::Vector{DateTime}
    anthropogenic::AnthropogenicInputs{FT}
    hvacschedule::HVACSchedule{FT}
    meteorological::MeteorologicalInputs{FT}
    sunposition::SunPositionInputs{FT}
end

function initialize_forcinginputset(
    ::Type{FT}, data::NCDataset, location::LocationProperties{FT}
) where {FT<:AbstractFloat}
    initialize(FT, ForcingInputSet, data, (FT,), location)
end

function TethysChlorisCore.get_required_fields(::Type{ForcingInputSet})
    return []
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{ForcingInputSet},
    data::NCDataset,
    params::Tuple,
    location::LocationProperties{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["datetime"] = Array(data["datetime"])

    processed["sunposition"] = initialize_sunposition_inputs(
        FT, data, processed["datetime"], location
    )
    processed["meteorological"] = initialize_meteorological_inputs(
        FT, data, processed["sunposition"].theta_Z
    )
    processed["anthropogenic"] = initialize_anthropogenic_inputs(
        FT, data, processed["meteorological"].Tatm
    )
    processed["hvacschedule"] = initialize_hvacschedule(FT, data)

    return processed
end

function TethysChlorisCore.validate_fields(::Type{ForcingInputSet}, data::NCDataset) end
