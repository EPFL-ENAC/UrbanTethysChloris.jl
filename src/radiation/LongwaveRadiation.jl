Base.@kwdef struct LongwaveRadiation{FT<:AbstractFloat}
    GroundImp::FT
    GroundBare::FT
    GroundVeg::FT
    Tree::FT
    WallSun::FT
    WallShade::FT
    TotalGround::FT
    TotalCanyon::FT
end

function Base.:-(
    a::LongwaveRadiation{FT}, b::LongwaveRadiation{FT}
) where {FT<:AbstractFloat}
    LongwaveRadiation{FT}(;
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

function combine(
    tree::LongwaveRadiation{FT}, notree::LongwaveRadiation{FT}, tree_fraction::FT
) where {FT<:AbstractFloat}
    LongwaveRadiation{FT}(;
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
