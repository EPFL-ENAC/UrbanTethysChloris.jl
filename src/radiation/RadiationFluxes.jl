abstract type AbstractRadiationFluxes{FT<:AbstractFloat} end

"""
    RadiationFluxes{FT<:AbstractFloat}

Structure representing radiation flux components in an urban canyon.

# Fields
- `GroundImp`: Radiation for impervious ground surface [W/m²]
- `GroundBare`: Radiation for bare ground surface [W/m²]
- `GroundVeg`: Radiation for vegetated ground surface [W/m²]
- `Tree`: Radiation for tree surface [W/m²]
- `WallSun`: Radiation for sunlit wall surface [W/m²]
- `WallShade`: Radiation for shaded wall surface [W/m²]
- `TotalGround`: Total radiation for all ground surfaces [W/m²]
- `TotalCanyon`: Total radiation for entire canyon [W/m²]
"""
Base.@kwdef struct RadiationFluxes{FT<:AbstractFloat} <: AbstractRadiationFluxes{FT}
    GroundImp::FT
    GroundBare::FT
    GroundVeg::FT
    Tree::FT
    WallSun::FT
    WallShade::FT
    TotalGround::FT
    TotalCanyon::FT
end

# TODO: replace by `diff`, instead of overloading `-`
function Base.:-(a::RadiationFluxes{FT}, b::RadiationFluxes{FT}) where {FT<:AbstractFloat}
    RadiationFluxes{FT}(;
        GroundImp=(a.GroundImp - b.GroundImp),
        GroundBare=(a.GroundBare - b.GroundBare),
        GroundVeg=(a.GroundVeg - b.GroundVeg),
        Tree=(a.Tree - b.Tree),
        WallSun=(a.WallSun - b.WallSun),
        WallShade=(a.WallShade - b.WallShade),
        TotalGround=(a.TotalGround - b.TotalGround),
        TotalCanyon=(a.TotalCanyon - b.TotalCanyon),
    )
end

function RadiationFluxes(::Type{FT}) where {FT<:AbstractFloat}
    RadiationFluxes{FT}(;
        GroundImp=zero(FT),
        GroundBare=zero(FT),
        GroundVeg=zero(FT),
        Tree=zero(FT),
        WallSun=zero(FT),
        WallShade=zero(FT),
        TotalGround=zero(FT),
        TotalCanyon=zero(FT),
    )
end

function RadiationFluxes(
    ::Type{FT}, data::AbstractDict, prefix::String=""
) where {FT<:AbstractFloat}
    return RadiationFluxes{FT}(;
        GroundImp=FT(data[prefix * "GroundImp"]),
        GroundBare=FT(data[prefix * "GroundBare"]),
        GroundVeg=FT(data[prefix * "GroundVeg"]),
        Tree=FT(data[prefix * "Tree"]),
        WallSun=FT(data[prefix * "WallSun"]),
        WallShade=FT(data[prefix * "WallShade"]),
        TotalGround=FT(data[prefix * "TotalGround"]),
        TotalCanyon=FT(data[prefix * "TotalCanyon"]),
    )
end

function Base.show(io::IO, obj::RadiationFluxes)
    print(io, typeof(obj))
    for field in fieldnames(typeof(obj))
        print(io, "\n", field, ": ", getfield(obj, field))
    end
end

Base.@kwdef struct AbsorbedRadiationFluxes{FT<:AbstractFloat} <: AbstractRadiationFluxes{FT}
    GroundImp::FT
    GroundBare::FT
    GroundVeg::FT
    Tree::FT
    WallSun::FT
    WallShade::FT
    TotalGround::FT
    TotalCanyon::FT
    absWindowsSun::FT
    transWindowSun::FT
    absWindowShade::FT
    transWindowShade::FT
    WallShadeExt::FT
    WallSunExt::FT
    WallShadeTransmitted::FT
    WallSunTransmitted::FT
end

function AbsorbedRadiationFluxes(
    x::RadiationFluxes{FT},
    absWindowsSun::FT,
    transWindowSun::FT,
    absWindowShade::FT,
    transWindowShade::FT,
    WallSunTransmitted::FT,
    WallShadeTransmitted::FT,
    WallSunExt::FT,
    WallShadeExt::FT,
) where {FT<:AbstractFloat}
    return AbsorbedRadiationFluxes{FT}(;
        GroundImp=x.GroundImp,
        GroundBare=x.GroundBare,
        GroundVeg=x.GroundVeg,
        Tree=x.Tree,
        WallSun=x.WallSun,
        WallShade=x.WallShade,
        TotalGround=x.TotalGround,
        TotalCanyon=x.TotalCanyon,
        absWindowsSun=absWindowsSun,
        transWindowSun=transWindowSun,
        absWindowShade=absWindowShade,
        transWindowShade=transWindowShade,
        WallShadeTransmitted=WallShadeTransmitted,
        WallSunTransmitted=WallSunTransmitted,
        WallSunExt=WallSunExt,
        WallShadeExt=WallShadeExt,
    )
end

"""
    interpolate(tree::RadiationFluxes{FT}, notree::RadiationFluxes{FT}, tree_fraction::FT) where {FT<:AbstractFloat}

Combines radiation components from scenarios with and without trees based on the tree fraction.

# Arguments
- `tree`: RadiationFluxes instance representing scenario with trees
- `notree`: RadiationFluxes instance representing scenario without trees
- `tree_fraction`: Fraction of area covered by trees [0-1]

# Returns
- `RadiationFluxes{FT}`: Combined radiation components
"""
function interpolate(
    tree::RadiationFluxes{FT}, notree::RadiationFluxes{FT}, tree_fraction::FT
) where {FT<:AbstractFloat}
    RadiationFluxes{FT}(;
        GroundImp=tree_fraction*tree.GroundImp + (1-tree_fraction)*notree.GroundImp,
        GroundBare=tree_fraction*tree.GroundBare + (1-tree_fraction)*notree.GroundBare,
        GroundVeg=tree_fraction*tree.GroundVeg + (1-tree_fraction)*notree.GroundVeg,
        Tree=tree.Tree,
        WallSun=tree_fraction*tree.WallSun + (1-tree_fraction)*notree.WallSun,
        WallShade=tree_fraction*tree.WallShade + (1-tree_fraction)*notree.WallShade,
        TotalGround=tree_fraction*tree.TotalGround + (1-tree_fraction)*notree.TotalGround,
        TotalCanyon=tree_fraction*tree.TotalCanyon + (1-tree_fraction)*notree.TotalCanyon,
    )
end
