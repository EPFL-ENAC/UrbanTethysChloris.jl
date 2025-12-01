"""
    ForcingInputSet{FT<:AbstractFloat, N} <: AbstractForcingInputSet{FT}

Forcing inputs for the Urban Tethys-Chloris model.

# Fields
- `datetime::Vector{DateTime}`: Date and time of the forcing inputs.
- `anthropogenic::AnthropogenicInputs{FT}`: Anthropogenic inputs.
- `hvacschedule::HVACSchedule{FT}`: HVAC schedule.
- `meteorological::MeteorologicalInputs{FT}`: Meteorological inputs.
- `sunposition::SunPositionInputs{FT}`: Sun position inputs.
"""
Base.@kwdef struct ForcingInputSet{FT<:AbstractFloat,N} <: Abstract1PForcingInputsSet{FT,N}
    datetime::Array{DateTime,N}
    anthropogenic::AnthropogenicInputs{FT,N}
    hvacschedule::HVACSchedule{FT,N}
    meteorological::MeteorologicalInputs{FT,N}
    sunposition::SunPositionInputs{FT,N}
end

function ForcingInputSet(
    ::Type{FT}, ::TimeSeries, data::NCDataset, location::LocationProperties{FT}
) where {FT<:AbstractFloat}
    initialize(FT, ForcingInputSet, data, (FT, dimension_value(TimeSeries())), location)
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

    processed["sunposition"] = SunPositionInputs(
        FT, TimeSeries(), data, processed["datetime"], location
    )
    processed["meteorological"] = MeteorologicalInputs(
        FT, TimeSeries(), data, processed["sunposition"].theta_Z
    )
    processed["anthropogenic"] = AnthropogenicInputs(
        FT, TimeSeries(), data, processed["meteorological"].Tatm
    )
    processed["hvacschedule"] = HVACSchedule(FT, TimeSeries(), data)

    return processed
end

function TethysChlorisCore.validate_fields(::Type{ForcingInputSet}, data::NCDataset) end

function Base.getindex(obj::T, idx::Int) where {FT<:AbstractFloat,T<:ForcingInputSet{FT,1}}
    scalar_type = typeof(obj).name.wrapper{FT,0}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        if field == :datetime
            fieldvals[field] = fill(getproperty(obj, field)[idx])
        else
            fieldvals[field] = getindex(getproperty(obj, field), idx)
        end
    end
    return scalar_type(; fieldvals...)
end
