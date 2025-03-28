abstract type AbstractParameterSet{FT<:AbstractFloat} <: AbstractModelComponentSet{FT} end
const APS = AbstractParameterSet

Base.@kwdef struct ParameterSet{FT<:AbstractFloat} <: APS{FT}
    soil::SoilParameters{FT}
    surfacefractions::SurfaceFractions{FT}
    thermal::ThermalProperties{FT}
    optical::OpticalProperties{FT}
    urbangeometry::UrbanGeometryParameters{FT}
    vegetation::VegetationParameters{FT}
end

function initialize_parameter_set(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, ParameterSet, data)
end

function TethysChlorisCore.get_required_fields(::Type{ParameterSet})
    return [:soil, :surfacefractions, :thermal, :optical, :urbangeometry, :vegetation]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{ParameterSet}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Initialize each component
    processed["soil"] = initialize_soil_parameters(FT, data["soil"])

    processed["surfacefractions"] = initialize_surfacefractions(
        FT, data["surfacefractions"]
    )
    processed["thermal"] = initialize_thermalproperties(FT, data["thermal"])
    processed["optical"] = initialize_optical_properties(
        FT, data["optical"], processed["surfacefractions"]
    )
    processed["urbangeometry"] = initialize_urbangeometry_parameters(
        FT, data["urbangeometry"]
    )
    processed["vegetation"] = initialize_vegetationparameters(FT, data["vegetation"])

    return processed
end
