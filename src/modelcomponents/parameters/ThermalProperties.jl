"""
    LocationSpecificThermalProperties{FT<:AbstractFloat} <: AbstractHeightDependentParameters{FT}

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
    TreeThermalProperties{FT<:AbstractFloat} <: AbstractHeightDependentParameters{FT}

Parameters for the tree thermal properties.

# Fields
- `Cthermal_leaf::FT`: [J m-2 K-1] Heat capacity per single leaf area
"""
Base.@kwdef struct TreeThermalProperties{FT<:AbstractFloat} <:
                   AbstractHeightDependentParameters{FT}
    Cthermal_leaf::FT
end

function initialize_tree_thermalproperties(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, TreeThermalProperties, data)
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
    tree::TreeThermalProperties{FT}
end

function initialize_thermalproperties(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, ThermalProperties, data)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{ThermalProperties}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    check_extraneous_fields(ThermalProperties, data)

    processed["roof"] = initialize(FT, LocationSpecificThermalProperties, data["roof"])
    processed["ground"] = initialize(FT, LocationSpecificThermalProperties, data["ground"])
    processed["wall"] = initialize(FT, LocationSpecificThermalProperties, data["wall"])
    processed["tree"] = initialize(FT, TreeThermalProperties, data["tree"])

    return processed
end
