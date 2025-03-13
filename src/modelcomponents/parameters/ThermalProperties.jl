"""
    LocationSpecificThermalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for the location-specific (wall, roof, ground) thermal properties.

# Fields
- `lan_dry::FT`: Thermal conductivity dry solid [W/m K]
- `cv_s::FT`: Volumetric heat capacity solid [J/m^3 K].
"""
Base.@kwdef struct LocationSpecificThermalProperties{FT<:AbstractFloat} <:
                   AbstractHeightDependentParameters{FT}
    lan_dry::FT
    cv_s::FT
end

function initialize_locationspecific_thermalproperties(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, LocationSpecificThermalProperties, data)
end

function TethysChlorisCore.get_required_fields(::Type{LocationSpecificThermalProperties})
    return [:lan_dry, :cv_s]
end

"""
    ThermalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}

Container for vegetation parameters for different urban surface components.

# Fields
- `roof::LocationSpecificThermalProperties{FT}`: Roof thermal properties
- `ground::LocationSpecificThermalProperties{FT}`: Ground thermal properties
- `wall::LocationSpecificThermalProperties{FT}`: Wall thermal properties
"""
Base.@kwdef struct ThermalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}
    roof::LocationSpecificThermalProperties{FT}
    ground::LocationSpecificThermalProperties{FT}
    wall::LocationSpecificThermalProperties{FT}
end

function initialize_thermalproperties(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, ThermalProperties, data)
end

function TethysChlorisCore.get_required_fields(::Type{ThermalProperties})
    return [:roof, :ground, :wall]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{ThermalProperties}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Check that data does not include a key beyond the three components
    for key in keys(data)
        if key ∉ ["roof", "ground", "wall"]
            throw(ArgumentError("Extraneous key: $key"))
        end
    end

    for (key, value) in data
        processed[key] = initialize(FT, LocationSpecificThermalProperties, value)
    end

    return processed
end
