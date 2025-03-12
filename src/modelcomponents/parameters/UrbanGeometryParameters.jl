"""
    UrbanGeometryParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for the UrbanGeometry model component.

# Fields
- `Height_canyon::FT`: Height of urban canyon [m].
- `Width_canyon::FT`: Ground width of urban canyon [m].
- `Width_roof::FT`: Roof width of urban canyon [m].
- `Height_tree::FT`: Tree height [m].
- `Radius_tree::FT`: Tree radius (=1/4 fg,tree *Wcan) [m].
- `Distance_tree::FT`: Distance of wall to tree trunk [m].
- `Hcan_max::FT`: Maximum height of roughness elements (buildings) [m].
- `Hcan_std::FT`: Standard deviation of roughness elements (buildings) [m].
- `trees::Bool`: Easy switch to include (=1) and exclude (=0) trees in the urban canyon.
- `ftree::FT`: Tree fraction along canyon axis
- `hcanyon::FT`: Normalized height of urban canyon [-].
- `wcanyon::FT`: Normalized ground width of urban canyon [-].
- `wroof::FT`: Normalized roof width of urban canyon [-]
- `htree::FT`: Normalized tree height [-]
- `radius_tree::FT`: Normalized tree radius [-].
- `distance_tree::FT`: Normalized distance of wall to tree trunk [-].
- `ratio::FT`: Height-to-width ratio [-].
- `wcanyon_norm::FT`: Normalized canyon width overall [-].
- `wroof_norm::FT`: Normalized roof width overall [-].
"""
Base.@kwdef struct UrbanGeometryParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    Height_canyon::FT
    Width_canyon::FT
    Width_roof::FT
    Height_tree::FT
    Radius_tree::FT
    Distance_tree::FT
    Hcan_max::FT
    Hcan_std::FT

    trees::Bool
    ftree::FT

    hcanyon::FT
    wcanyon::FT
    wroof::FT
    htree::FT
    radius_tree::FT
    distance_tree::FT
    ratio::FT

    wcanyon_norm::FT
    wroof_norm::FT
end

function initialize_urbangeometry_parameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, UrbanGeometryParameters, data)
end

function TethysChlorisCore.get_required_fields(::Type{UrbanGeometryParameters})
    return [
        :Height_canyon,
        :Width_canyon,
        :Width_roof,
        :Height_tree,
        :Radius_tree,
        :Distance_tree,
        :Hcan_max,
        :Hcan_std,
        :trees,
        :ftree,
    ]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{UrbanGeometryParameters}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    data["hcanyon"] = data["Height_canyon"] / data["Width_canyon"]
    data["wcanyon"] = data["Width_canyon"] / data["Width_canyon"]
    data["wroof"] = data["Width_roof"] / data["Width_canyon"]
    data["htree"] = data["Height_tree"] / data["Width_canyon"]
    data["radius_tree"] = data["Radius_tree"] / data["Width_canyon"]
    data["distance_tree"] = data["Distance_tree"] / data["Width_canyon"]
    data["ratio"] = data["hcanyon"] / data["wcanyon"]
    data["wcanyon_norm"] = data["wcanyon"] / (data["wcanyon"] + data["wroof"])
    data["wroof_norm"] = data["wroof"] / (data["wcanyon"] + data["wroof"])

    return data
end
