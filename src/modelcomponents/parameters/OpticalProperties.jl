"""
    SimpleOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for the location-specific (wall and tree) optical properties, for which albedo
and emissivity have been pre-calculated.

# Fields
- `albedo::FT`: Surface albedo (-)
- `emissivity::FT`: Surface emissivity (-)
"""
Base.@kwdef struct SimpleOpticalProperties{FT<:AbstractFloat} <:
                   AbstractHeightDependentParameters{FT}
    albedo::FT
    emissivity::FT
end

function initialize_simple_opticalproperties(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, SimpleOpticalProperties, data)
end

function SimpleOpticalProperties(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return SimpleOpticalProperties{FT}(data["albedo"], data["emissivity"])
end

"""
    VegetatedOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}

Optical properties specific to vegetated surfaces (roof and ground)

# Fields
- `aveg::FT`: Vegetation surface albedo (-)
- `aimp::FT`: Impervious surface albedo (-)
- `abare::FT`: Bare surface albedo (-), only used for ground
- `eveg::FT`: Vegetation surface emissivity (-)
- `eimp::FT`: Impervious surface emissivity (-)
- `ebare::FT`: Bare surface emissivity (-), only used for ground
"""
Base.@kwdef struct VegetatedOpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}
    aveg::FT
    aimp::FT
    abare::FT = FT(NaN)  # Only used for ground
    albedo::FT

    eveg::FT
    eimp::FT
    ebare::FT = FT(NaN)  # Only used for ground
    emissivity::FT
end

function initialize_vegetated_opticalproperties(
    ::Type{FT}, data::Dict{String,Any}, fractions::LocationSpecificSurfaceFractions{FT}
) where {FT<:AbstractFloat}
    return initialize(FT, VegetatedOpticalProperties, data, (FT,), fractions)
end

function VegetatedOpticalProperties(
    ::Type{FT}, data::AbstractDict
) where {FT<:AbstractFloat}
    return VegetatedOpticalProperties{FT}(
        data["aveg"],
        data["aimp"],
        get(data, "abare", FT(NaN)),
        data["albedo"],
        data["eveg"],
        data["eimp"],
        get(data, "ebare", FT(NaN)),
        data["emissivity"],
    )
end

function TethysChlorisCore.get_optional_fields(::Type{VegetatedOpticalProperties})
    return [:abare, :ebare]
end

function TethysChlorisCore.get_calculated_fields(::Type{VegetatedOpticalProperties})
    return [:albedo, :emissivity]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{VegetatedOpticalProperties},
    data::Dict{String,Any},
    params::Tuple,
    fractions::LocationSpecificSurfaceFractions{FT},
) where {FT<:AbstractFloat}
    processed = copy(data)

    processed["albedo"] = fractions.fveg * data["aveg"] + fractions.fimp * data["aimp"]

    if haskey(data, "abare") && !isnan(data["abare"])
        processed["albedo"] += fractions.fbare * data["abare"]
    end

    processed["emissivity"] = fractions.fveg * data["eveg"] + fractions.fimp * data["eimp"]

    if haskey(data, "ebare") && !isnan(data["ebare"])
        processed["emissivity"] += fractions.fbare * data["ebare"]
    end

    return processed
end

"""
    OpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}

Container for optical properties for different urban surface components.

# Fields
- `roof::VegetatedOpticalProperties{FT}`: Roof optical properties
- `ground::VegetatedOpticalProperties{FT}`: Ground optical properties
- `wall::SimpleOpticalProperties{FT}`: Wall optical properties
- `tree::SimpleOpticalProperties{FT}`: Tree optical properties
"""
Base.@kwdef struct OpticalProperties{FT<:AbstractFloat} <: AbstractParameters{FT}
    roof::VegetatedOpticalProperties{FT}
    ground::VegetatedOpticalProperties{FT}
    wall::SimpleOpticalProperties{FT}
    tree::SimpleOpticalProperties{FT}
end

function initialize_optical_properties(
    ::Type{FT}, data::Dict{String,Any}, fractions::SurfaceFractions{FT}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["roof"] = initialize_vegetated_opticalproperties(
        FT, data["roof"], fractions.roof
    )
    processed["ground"] = initialize_vegetated_opticalproperties(
        FT, data["ground"], fractions.ground
    )
    processed["wall"] = initialize_simple_opticalproperties(FT, data["wall"])
    processed["tree"] = initialize_simple_opticalproperties(FT, data["tree"])

    return initialize(FT, OpticalProperties, processed)
end
