abstract type AbstractExtrapoledModelVariables{FT<:AbstractFloat} end

"""
    extrapolate!(y::AbstractExtrapoledModelVariables{FT}, x::T, i::Signed) where {FT<:AbstractFloat,T}

Extrapolates the values in `x` into the extrapolated model variables `y`. The parameter `i`
indicates the current iteration number. If `i > 2`, a linear extrapolation is performed
based on the last two known values. If `i ≤ 2`, the values from `x` are simply duplicated to
fill the extrapolated structure.
"""
function extrapolate!(
    y::AbstractExtrapoledModelVariables{FT}, x::T, i::Signed
) where {FT<:AbstractFloat,T}
    if i > 2
        extrapolate!(y, x)
    else
        for field in fieldnames(T)
            tm1 = getfield(x, field)
            setfield!(y, field, SVector{4,FT}(tm1, tm1, tm1, tm1))
        end
    end
end

function extrapolate!(
    y::AbstractExtrapoledModelVariables{FT}, x::T
) where {FT<:AbstractFloat,T}
    for field in fieldnames(T)
        tm1 = getfield(x, field)
        tm2 = getfield(y, field)[2]
        diff = tm1 - tm2

        setfield!(y, field, SVector{4,FT}(tm2, tm1, tm1 + diff, tm1 + 2*diff))
    end
end

function Base.show(io::IO, obj::AbstractExtrapoledModelVariables)
    print(io, typeof(obj))
    for field in fieldnames(typeof(obj))
        print(io, "\n", field, ": ", getfield(obj, field))
    end
end

Base.@kwdef mutable struct ExtrapolatedTempVec{FT<:AbstractFloat} <:
                           AbstractExtrapoledModelVariables{FT}
    TRoofImp::SVector{4,FT}
    TRoofVeg::SVector{4,FT}
    TRoofIntImp::SVector{4,FT}
    TRoofIntVeg::SVector{4,FT}
    TGroundImp::SVector{4,FT}
    TGroundBare::SVector{4,FT}
    TGroundVeg::SVector{4,FT}
    TTree::SVector{4,FT}
    TWallSun::SVector{4,FT}
    TWallShade::SVector{4,FT}
    TWallIntSun::SVector{4,FT}
    TWallIntShade::SVector{4,FT}
    TCanyon::SVector{4,FT}
    Tatm::SVector{4,FT}
    T2m::SVector{4,FT}
end

function ExtrapolatedTempVec(::Type{FT}) where {FT<:AbstractFloat}
    return ExtrapolatedTempVec{FT}(;
        TRoofImp=zero(SVector{4,FT}),
        TRoofVeg=zero(SVector{4,FT}),
        TRoofIntImp=zero(SVector{4,FT}),
        TRoofIntVeg=zero(SVector{4,FT}),
        TGroundImp=zero(SVector{4,FT}),
        TGroundBare=zero(SVector{4,FT}),
        TGroundVeg=zero(SVector{4,FT}),
        TTree=zero(SVector{4,FT}),
        TWallSun=zero(SVector{4,FT}),
        TWallShade=zero(SVector{4,FT}),
        TWallIntSun=zero(SVector{4,FT}),
        TWallIntShade=zero(SVector{4,FT}),
        TCanyon=zero(SVector{4,FT}),
        Tatm=zero(SVector{4,FT}),
        T2m=zero(SVector{4,FT}),
    )
end

function ExtrapolatedTempVec(
    x::ModelComponents.ModelVariables.TempVec{FT}
) where {FT<:AbstractFloat}
    return ExtrapolatedTempVec{FT}(;
        TRoofImp=SVector{4,FT}(x.TRoofImp, x.TRoofImp, x.TRoofImp, x.TRoofImp),
        TRoofVeg=SVector{4,FT}(x.TRoofVeg, x.TRoofVeg, x.TRoofVeg, x.TRoofVeg),
        TRoofIntImp=SVector{4,FT}(
            x.TRoofIntImp, x.TRoofIntImp, x.TRoofIntImp, x.TRoofIntImp
        ),
        TRoofIntVeg=SVector{4,FT}(
            x.TRoofIntVeg, x.TRoofIntVeg, x.TRoofIntVeg, x.TRoofIntVeg
        ),
        TGroundImp=SVector{4,FT}(x.TGroundImp, x.TGroundImp, x.TGroundImp, x.TGroundImp),
        TGroundBare=SVector{4,FT}(
            x.TGroundBare, x.TGroundBare, x.TGroundBare, x.TGroundBare
        ),
        TGroundVeg=SVector{4,FT}(x.TGroundVeg, x.TGroundVeg, x.TGroundVeg, x.TGroundVeg),
        TTree=SVector{4,FT}(x.TTree, x.TTree, x.TTree, x.TTree),
        TWallSun=SVector{4,FT}(x.TWallSun, x.TWallSun, x.TWallSun, x.TWallSun),
        TWallShade=SVector{4,FT}(x.TWallShade, x.TWallShade, x.TWallShade, x.TWallShade),
        TWallIntSun=SVector{4,FT}(
            x.TWallIntSun, x.TWallIntSun, x.TWallIntSun, x.TWallIntSun
        ),
        TWallIntShade=SVector{4,FT}(
            x.TWallIntShade, x.TWallIntShade, x.TWallIntShade, x.TWallIntShade
        ),
        TCanyon=SVector{4,FT}(x.TCanyon, x.TCanyon, x.TCanyon, x.TCanyon),
        Tatm=SVector{4,FT}(x.Tatm, x.Tatm, x.Tatm, x.Tatm),
        T2m=SVector{4,FT}(x.T2m, x.T2m, x.T2m, x.T2m),
    )
end

Base.@kwdef mutable struct ExtrapolatedHumidity{FT} <: AbstractExtrapoledModelVariables{FT}
    CanyonRelative::SVector{4,FT}
    CanyonSpecific::SVector{4,FT}
    CanyonVapourPre::SVector{4,FT}
    CanyonRelativeSat::SVector{4,FT}
    CanyonSpecificSat::SVector{4,FT}
    CanyonVapourPreSat::SVector{4,FT}
    AtmRelative::SVector{4,FT}
    AtmSpecific::SVector{4,FT}
    AtmVapourPre::SVector{4,FT}
    AtmRelativeSat::SVector{4,FT}
    AtmSpecificSat::SVector{4,FT}
    AtmVapourPreSat::SVector{4,FT}
    q2m::SVector{4,FT}
end

function ExtrapolatedHumidity(::Type{FT}) where {FT<:AbstractFloat}
    return ExtrapolatedHumidity{FT}(;
        CanyonRelative=zero(SVector{4,FT}),
        CanyonSpecific=zero(SVector{4,FT}),
        CanyonVapourPre=zero(SVector{4,FT}),
        CanyonRelativeSat=zero(SVector{4,FT}),
        CanyonSpecificSat=zero(SVector{4,FT}),
        CanyonVapourPreSat=zero(SVector{4,FT}),
        AtmRelative=zero(SVector{4,FT}),
        AtmSpecific=zero(SVector{4,FT}),
        AtmVapourPre=zero(SVector{4,FT}),
        AtmRelativeSat=zero(SVector{4,FT}),
        AtmSpecificSat=zero(SVector{4,FT}),
        AtmVapourPreSat=zero(SVector{4,FT}),
        q2m=zero(SVector{4,FT}),
    )
end

function ExtrapolatedHumidity(
    x::ModelComponents.ModelVariables.Humidity{FT}
) where {FT<:AbstractFloat}
    return ExtrapolatedHumidity{FT}(;
        CanyonRelative=SVector{4,FT}(
            x.CanyonRelative, x.CanyonRelative, x.CanyonRelative, x.CanyonRelative
        ),
        CanyonSpecific=SVector{4,FT}(
            x.CanyonSpecific, x.CanyonSpecific, x.CanyonSpecific, x.CanyonSpecific
        ),
        CanyonVapourPre=SVector{4,FT}(
            x.CanyonVapourPre, x.CanyonVapourPre, x.CanyonVapourPre, x.CanyonVapourPre
        ),
        CanyonRelativeSat=SVector{4,FT}(
            x.CanyonRelativeSat,
            x.CanyonRelativeSat,
            x.CanyonRelativeSat,
            x.CanyonRelativeSat,
        ),
        CanyonSpecificSat=SVector{4,FT}(
            x.CanyonSpecificSat,
            x.CanyonSpecificSat,
            x.CanyonSpecificSat,
            x.CanyonSpecificSat,
        ),
        CanyonVapourPreSat=SVector{4,FT}(
            x.CanyonVapourPreSat,
            x.CanyonVapourPreSat,
            x.CanyonVapourPreSat,
            x.CanyonVapourPreSat,
        ),
        AtmRelative=SVector{4,FT}(
            x.AtmRelative, x.AtmRelative, x.AtmRelative, x.AtmRelative
        ),
        AtmSpecific=SVector{4,FT}(
            x.AtmSpecific, x.AtmSpecific, x.AtmSpecific, x.AtmSpecific
        ),
        AtmVapourPre=SVector{4,FT}(
            x.AtmVapourPre, x.AtmVapourPre, x.AtmVapourPre, x.AtmVapourPre
        ),
        AtmRelativeSat=SVector{4,FT}(
            x.AtmRelativeSat, x.AtmRelativeSat, x.AtmRelativeSat, x.AtmRelativeSat
        ),
        AtmSpecificSat=SVector{4,FT}(
            x.AtmSpecificSat, x.AtmSpecificSat, x.AtmSpecificSat, x.AtmSpecificSat
        ),
        AtmVapourPreSat=SVector{4,FT}(
            x.AtmVapourPreSat, x.AtmVapourPreSat, x.AtmVapourPreSat, x.AtmVapourPreSat
        ),
        q2m=SVector{4,FT}(x.q2m, x.q2m, x.q2m, x.q2m),
    )
end

Base.@kwdef mutable struct ExtrapolatedTempVecB{FT<:AbstractFloat} <:
                           AbstractExtrapoledModelVariables{FT}
    Tceiling::SVector{4,FT}
    Tinwallsun::SVector{4,FT}
    Tinwallshd::SVector{4,FT}
    Twindows::SVector{4,FT}
    Tinground::SVector{4,FT}
    Tintmass::SVector{4,FT}
    Tbin::SVector{4,FT}
    qbin::SVector{4,FT}
end

function ExtrapolatedTempVecB(::Type{FT}) where {FT<:AbstractFloat}
    return ExtrapolatedTempVecB{FT}(;
        Tceiling=zero(SVector{4,FT}),
        Tinwallsun=zero(SVector{4,FT}),
        Tinwallshd=zero(SVector{4,FT}),
        Twindows=zero(SVector{4,FT}),
        Tinground=zero(SVector{4,FT}),
        Tintmass=zero(SVector{4,FT}),
        Tbin=zero(SVector{4,FT}),
        qbin=zero(SVector{4,FT}),
    )
end

function ExtrapolatedTempVecB(
    x::ModelComponents.ModelVariables.TempVecB{FT}
) where {FT<:AbstractFloat}
    return ExtrapolatedTempVecB{FT}(;
        Tceiling=SVector{4,FT}(x.Tceiling, x.Tceiling, x.Tceiling, x.Tceiling),
        Tinwallsun=SVector{4,FT}(x.Tinwallsun, x.Tinwallsun, x.Tinwallsun, x.Tinwallsun),
        Tinwallshd=SVector{4,FT}(x.Tinwallshd, x.Tinwallshd, x.Tinwallshd, x.Tinwallshd),
        Twindows=SVector{4,FT}(x.Twindows, x.Twindows, x.Twindows, x.Twindows),
        Tinground=SVector{4,FT}(x.Tinground, x.Tinground, x.Tinground, x.Tinground),
        Tintmass=SVector{4,FT}(x.Tintmass, x.Tintmass, x.Tintmass, x.Tintmass),
        Tbin=SVector{4,FT}(x.Tbin, x.Tbin, x.Tbin, x.Tbin),
        qbin=SVector{4,FT}(x.qbin, x.qbin, x.qbin, x.qbin),
    )
end

Base.@kwdef mutable struct Meteotm1{FT<:AbstractFloat}
    SWRin::SVector{2,FT}
    Rain::SVector{2,FT}
end

function Meteotm1(
    x::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0}
) where {FT<:AbstractFloat}
    SWRin = x.SAB1_in + x.SAB2_in + x.SAD1_in + x.SAD2_in
    return Meteotm1{FT}(;
        SWRin=SVector{2,FT}(SWRin, SWRin), Rain=SVector{2,FT}(x.Rain, x.Rain)
    )
end

function update!(
    y::Meteotm1{FT}, x::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0}
) where {FT<:AbstractFloat}
    SWRin = x.SAB1_in + x.SAB2_in + x.SAD1_in + x.SAD2_in
    y.SWRin = SVector{2,FT}(y.SWRin[2], SWRin)
    y.Rain = SVector{2,FT}(y.Rain[2], x.Rain)

    return nothing
end
