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
