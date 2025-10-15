"""
    HVACSchedule{FT<:AbstractFloat}

HVAC schedule parameters that specify heat and humidity sources from building equipment and occupants.

# Fields
- `Hequip`: Sensible heat from equipment [W/m² building ground area]
- `Hpeople`: Sensible heat from people [W/m² building ground area]
- `LEequip`: Latent heat from equipment [W/m² building ground area]
- `LEpeople`: Latent heat from people [W/m² building ground area]
- `AirConRoomFraction`: Fraction of air conditioned space [-]
"""
Base.@kwdef struct HVACSchedule{FT<:AbstractFloat} <: AbstractParameters{FT}
    Hequip::Vector{FT}
    Hpeople::Vector{FT}
    LEequip::Vector{FT}
    LEpeople::Vector{FT}
    AirConRoomFraction::Vector{FT}
end

function TethysChlorisCore.get_optional_fields(::Type{HVACSchedule})
    return [:Hequip, :Hpeople, :LEequip, :LEpeople, :AirConRoomFraction]
end

function initialize_hvacschedule(::Type{FT}, data::NCDataset) where {FT<:AbstractFloat}
    return initialize(FT, HVACSchedule, data, (FT,))
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{HVACSchedule}, data::NCDataset, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    for field in TethysChlorisCore.get_optional_fields(HVACSchedule)
        if haskey(data, field)
            if dimnames(data[field]) == ()
                processed[String(field)] = fill(FT(data[field][]), data.dim["hours"])
            else
                processed[String(field)] = Array(data[field])
            end
        else
            if field == :AirConRoomFraction
                processed[String(field)] = ones(FT, data.dim["hours"])
            else
                processed[String(field)] = zeros(FT, data.dim["hours"])
            end
        end
    end

    return processed
end

function TethysChlorisCore.validate_fields(::Type{HVACSchedule}, data::NCDataset) end
