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
Base.@kwdef struct RadiationFluxes{FT<:AbstractFloat}
    GroundImp::FT
    GroundBare::FT
    GroundVeg::FT
    Tree::FT
    WallSun::FT
    WallShade::FT
    TotalGround::FT
    TotalCanyon::FT
end

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
