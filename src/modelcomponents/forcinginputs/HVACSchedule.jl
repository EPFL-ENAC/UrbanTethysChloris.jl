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
    Hequip::FT
    Hpeople::FT
    LEequip::FT
    LEpeople::FT
    AirConRoomFraction::FT
end

"""
    initialize_hvacschedule(::Type{FT}, data::Dict{String,Any}) where {FT<:AbstractFloat}

Initialize a `HVACSchedule` instance from a dictionary of parameters.

# Arguments
- `FT`: Float type
- `data`: Dictionary containing HVAC schedule parameters
"""
function initialize_hvacschedule(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, HVACSchedule, data)
end
