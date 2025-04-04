Base.@kwdef struct LocationSpecificSurfaceFractions{FT<:AbstractFloat} <:
                   AbstractHeightDependentParameters{FT}
    fveg::FT
    fbare::FT = FT(NaN)
    fimp::FT
    Per_runoff::FT
end

function initialize_locationspecific_surfacefractions(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, LocationSpecificSurfaceFractions, data)
end

function get_optional_fields(::Type{LocationSpecificSurfaceFractions})
    return [:fbare]
end

function TethysChlorisCore.validate_fields(
    ::Type{LocationSpecificSurfaceFractions}, data::Dict{String,Any}
)
    check_extraneous_fields(LocationSpecificSurfaceFractions, data)

    if data["fveg"] + data["fimp"] != 1.0
        throw(ArgumentError("Surface fractions must sum to 1.0"))
    end

    # if bare is a field of data and not NaN
    if haskey(data, "fbare") && !isnan(data["fbare"])
        if data["fveg"] + data["fbare"] + data["fimp"] != 1.0
            throw(ArgumentError("Surface fractions must sum to 1.0"))
        end
    end

    if data["Per_runoff"] < 0.0 || data["Per_runoff"] > 1.0
        throw(ArgumentError("Runoff fraction must be between 0.0 and 1.0"))
    end
end

"""
    SurfaceFractions{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for the SurfaceFractions model component.

# Fields
- `roof::LocationSpecificSurfaceFractions{FT}`: Roof surface fractions.
- `ground::LocationSpecificSurfaceFractions{FT}`: Ground surface fractions.
"""
Base.@kwdef struct SurfaceFractions{FT<:AbstractFloat} <: AbstractParameters{FT}
    roof::LocationSpecificSurfaceFractions{FT}
    ground::LocationSpecificSurfaceFractions{FT}
end

function initialize_surfacefractions(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, SurfaceFractions, data)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{SurfaceFractions}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["roof"] = initialize_locationspecific_surfacefractions(FT, data["roof"])
    processed["ground"] = initialize_locationspecific_surfacefractions(FT, data["ground"])

    return processed
end
