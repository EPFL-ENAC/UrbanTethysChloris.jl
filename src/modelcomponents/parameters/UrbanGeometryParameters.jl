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
    return initialize(FT, UrbanGeometryParameters, data, (FT,))
end

function UrbanGeometryParameters(FT, Gemeotry_m::AbstractDict)
    return UrbanGeometryParameters{FT}(;
        Height_canyon=FT(Gemeotry_m["Height_canyon"]),
        Width_canyon=FT(Gemeotry_m["Width_canyon"]),
        Width_roof=FT(Gemeotry_m["Width_roof"]),
        Height_tree=FT(Gemeotry_m["Height_tree"]),
        Radius_tree=FT(Gemeotry_m["Radius_tree"]),
        Distance_tree=FT(Gemeotry_m["Distance_tree"]),
        Hcan_max=FT(Gemeotry_m["Hcan_max"]),
        Hcan_std=FT(Gemeotry_m["Hcan_std"]),
        trees=false,
        ftree=FT(NaN),
        hcanyon=FT(NaN),
        wcanyon=FT(NaN),
        wroof=FT(NaN),
        htree=FT(NaN),
        radius_tree=FT(NaN),
        distance_tree=FT(NaN),
        ratio=FT(NaN),
        wcanyon_norm=FT(NaN),
        wroof_norm=FT(NaN),
    )
end

function UrbanGeometryParameters(
    FT, Gemeotry_m::AbstractDict, geometry::AbstractDict, ParTree::AbstractDict
)
    return UrbanGeometryParameters{FT}(;
        Height_canyon=FT(Gemeotry_m["Height_canyon"]),
        Width_canyon=FT(Gemeotry_m["Width_canyon"]),
        Width_roof=FT(Gemeotry_m["Width_roof"]),
        Height_tree=FT(Gemeotry_m["Height_tree"]),
        Radius_tree=FT(Gemeotry_m["Radius_tree"]),
        Distance_tree=FT(Gemeotry_m["Distance_tree"]),
        Hcan_max=FT(Gemeotry_m["Hcan_max"]),
        Hcan_std=FT(Gemeotry_m["Hcan_std"]),
        trees=Bool(ParTree["trees"]),
        ftree=FT(ParTree["ftree"]),
        hcanyon=FT(geometry["hcanyon"]),
        wcanyon=FT(geometry["wcanyon"]),
        wroof=FT(geometry["wroof"]),
        htree=FT(geometry["htree"]),
        radius_tree=FT(geometry["radius_tree"]),
        distance_tree=FT(geometry["distance_tree"]),
        ratio=FT(geometry["ratio"]),
        wcanyon_norm=FT(geometry["wcanyon_norm"]),
        wroof_norm=FT(geometry["wroof_norm"]),
    )
end

function TethysChlorisCore.get_calculated_fields(::Type{UrbanGeometryParameters})
    return [
        :hcanyon,
        :wcanyon,
        :wroof,
        :htree,
        :radius_tree,
        :distance_tree,
        :ratio,
        :wcanyon_norm,
        :wroof_norm,
    ]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{UrbanGeometryParameters}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = copy(data)

    hcanyon, wcanyon, wroof, htree, radius_tree, distance_tree, ratio, wcanyon_norm, wroof_norm = preprocess_geometry(
        data["Height_canyon"],
        data["Width_canyon"],
        data["Width_roof"],
        data["Height_tree"],
        data["Radius_tree"],
        data["Distance_tree"],
    )

    processed["hcanyon"] = hcanyon
    processed["wcanyon"] = wcanyon
    processed["wroof"] = wroof
    processed["htree"] = htree
    processed["radius_tree"] = radius_tree
    processed["distance_tree"] = distance_tree
    processed["ratio"] = ratio
    processed["wcanyon_norm"] = wcanyon_norm
    processed["wroof_norm"] = wroof_norm

    return processed
end

function TethysChlorisCore.validate_fields(
    ::Type{UrbanGeometryParameters}, data::Dict{String,Any}
)
    check_extraneous_fields(UrbanGeometryParameters, data)

    # check that none of the data fields used by the validate function is NaN
    for key in keys(data)
        if key in [
            "Height_canyon",
            "Width_canyon",
            "Width_roof",
            "Height_tree",
            "Radius_tree",
            "Distance_tree",
        ]
            if isnan(data[key])
                throw(ArgumentError("Field $key must be a number, got NaN"))
            end
        end
    end

    if data["Width_canyon"] <= 0.0
        throw(ArgumentError("Urban canyon width must be strictly positive"))
    end

    if data["Width_roof"] <= 0.0
        throw(ArgumentError("Urban canyon roof width must be strictly positive"))
    end

    if data["Height_canyon"] <= 0.0
        throw(ArgumentError("Urban canyon height must be strictly positive"))
    end

    if data["ftree"] != 1.0
        throw(ArgumentError("Tree fraction must be 1.0"))
    end

    hcanyon, wcanyon, wroof, htree, radius_tree, distance_tree = preprocess_geometry(
        data["Height_canyon"],
        data["Width_canyon"],
        data["Width_roof"],
        data["Height_tree"],
        data["Radius_tree"],
        data["Distance_tree"],
    )

    if data["trees"]
        if data["Height_tree"] <= 0.0
            throw(ArgumentError("Tree height cannot be 0"))
        end

        if data["Radius_tree"] <= 0.0
            throw(ArgumentError("Tree radius cannot be 0"))
        end

        if data["Distance_tree"] <= 0.0
            throw(ArgumentError("Tree distance cannot be 0"))
        end

        if hcanyon - 2 * radius_tree <= 0.0
            throw(ArgumentError("Tree radius is too large for the canyon height"))
        end

        if wcanyon - 4 * radius_tree <= 0.0
            throw(ArgumentError("Tree radius is too large for the canyon width"))
        end

        if hcanyon - (radius_tree + htree) <= 0.0
            throw(ArgumentError("Tree height is too large for the canyon height"))
        end

        if htree - radius_tree <= 0.0
            throw(ArgumentError("Tree radius is too large for the tree height"))
        end

        if distance_tree - radius_tree <= 0.0
            throw(ArgumentError("Tree radius is too large for the tree distance"))
        end

        if wcanyon - 2 * (radius_tree + distance_tree) <= 0.0
            throw(ArgumentError("Tree distance is too large for the canyon width"))
        end
    end

    if hcanyon <= 0.0
        throw(ArgumentError("Normalized height of urban canyon must be strictly positive"))
    end

    if wcanyon < 0.0
        throw(ArgumentError("Normalized width of urban canyon must be stricly positive"))
    end

    if wroof < 0.0
        throw(
            ArgumentError("Normalized roof width of urban canyon must be strictly positive")
        )
    end
end

function preprocess_geometry(
    canyon_height::FT,
    canon_width::FT,
    roof_width::FT,
    tree_height::FT,
    tree_radius::FT,
    tree_distance::FT,
) where {FT<:AbstractFloat}
    hcanyon = canyon_height / canon_width
    wcanyon = canon_width / canon_width
    wroof = roof_width / canon_width
    htree = tree_height / canon_width
    radius_tree = tree_radius / canon_width
    distance_tree = tree_distance / canon_width
    ratio = hcanyon / wcanyon
    wcanyon_norm = wcanyon / (wcanyon + wroof)
    wroof_norm = wroof / (wcanyon + wroof)

    return hcanyon,
    wcanyon, wroof, htree, radius_tree, distance_tree, ratio, wcanyon_norm,
    wroof_norm
end
