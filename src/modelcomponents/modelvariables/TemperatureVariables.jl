"""
    TempVec{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Temperature Vector Fields.

# Fields
- `TRoofImp`: Temperature roof impervious area [K]
- `TRoofVeg`: Temperature roof vegetated area [K]
- `TRoofIntImp`: Interior temperature roof impervious area [K]
- `TRoofIntVeg`: Interior temperature roof vegetated area [K]
- `TGroundImp`: Temperature ground impervious area [K]
- `TGroundBare`: Temperature ground bare area [K]
- `TGroundVeg`: Temperature ground vegetated area [K]
- `TTree`: Temperature tree canopy [K]
- `TWallSun`: Temperature sunlit wall area [K]
- `TWallShade`: Temperature shaded wall area [K]
- `TWallIntSun`: Interior temperature sunlit wall [K]
- `TWallIntShade`: Interior temperature shaded wall [K]
- `TCanyon`: Temperature canyon [K]
- `Tatm`: Temperature atmosphere (measured) [K]
"""
Base.@kwdef mutable struct TempVec{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    TRoofImp::FT
    TRoofVeg::FT
    TRoofIntImp::FT
    TRoofIntVeg::FT
    TGroundImp::FT
    TGroundBare::FT
    TGroundVeg::FT
    TTree::FT
    TWallSun::FT
    TWallShade::FT
    TWallIntSun::FT
    TWallIntShade::FT
    TCanyon::FT
    Tatm::FT
    T2m::FT
end

function TempVec(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, TempVec, Dict{String,Any}())
end

function TempVec(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return TempVec{FT}(;
        TRoofImp=data["TRoofImp"],
        TRoofVeg=data["TRoofVeg"],
        TRoofIntImp=data["TRoofIntImp"],
        TRoofIntVeg=data["TRoofIntVeg"],
        TGroundImp=data["TGroundImp"],
        TGroundBare=data["TGroundBare"],
        TGroundVeg=data["TGroundVeg"],
        TTree=data["TTree"],
        TWallSun=data["TWallSun"],
        TWallShade=data["TWallShade"],
        TWallIntSun=data["TWallIntSun"],
        TWallIntShade=data["TWallIntShade"],
        TCanyon=data["TCanyon"],
        Tatm=data["Tatm"],
        T2m=data["T2m"],
    )
end

function update!(x::TempVec{FT}, Ttot::Vector{FT}) where {FT<:AbstractFloat}
    x.TRoofImp = Ttot[1]
    x.TRoofVeg = Ttot[2]
    x.TRoofIntImp = Ttot[3]
    x.TRoofIntVeg = Ttot[4]
    x.TGroundImp = Ttot[5]
    x.TGroundBare = Ttot[6]
    x.TGroundVeg = Ttot[7]
    x.TWallSun = Ttot[8]
    x.TWallShade = Ttot[9]
    x.TTree = Ttot[10]
    x.TWallIntSun = Ttot[11]
    x.TWallIntShade = Ttot[12]
    x.TCanyon = Ttot[13]

    return nothing
end

function update!(x::TempVec{FT}, y::TempVec{FT}) where {FT<:AbstractFloat}
    x.TRoofImp = y.TRoofImp
    x.TRoofVeg = y.TRoofVeg
    x.TRoofIntImp = y.TRoofIntImp
    x.TRoofIntVeg = y.TRoofIntVeg
    x.TGroundImp = y.TGroundImp
    x.TGroundBare = y.TGroundBare
    x.TGroundVeg = y.TGroundVeg
    x.TTree = y.TTree
    x.TWallSun = y.TWallSun
    x.TWallShade = y.TWallShade
    x.TWallIntSun = y.TWallIntSun
    x.TWallIntShade = y.TWallIntShade
    x.TCanyon = y.TCanyon
    x.Tatm = y.Tatm
    x.T2m = y.T2m
end

"""
    TempDamp{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Temperature Dampening Fields.

# Fields
- `TDampGroundImp`: Dampening temperature ground impervious area [K]
- `TDampGroundBare`: Dampening temperature ground bare area [K]
- `TDampGroundVeg`: Dampening temperature ground vegetated area [K]
- `TDampTree`: Dampening temperature tree canopy [K]
- `TDampGroundBuild`: Dampening temperature of ground in building interior [K]
"""
Base.@kwdef mutable struct TempDamp{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    TDampGroundImp::FT
    TDampGroundBare::FT
    TDampGroundVeg::FT
    TDampTree::FT
    TDampGroundBuild::FT
end

function TempDamp(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, TempDamp, Dict{String,Any}())
end

function TempDamp(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return TempDamp{FT}(;
        TDampGroundImp=data["TDampGroundImp"],
        TDampGroundBare=data["TDampGroundBare"],
        TDampGroundVeg=data["TDampGroundVeg"],
        TDampTree=data["TDampTree"],
        TDampGroundBuild=data["TDampGroundBuild"],
    )
end

# Necessary to avoid a reference assignment
function update!(x::TempDamp{FT}, y::TempDamp{FT}) where {FT<:AbstractFloat}
    x.TDampGroundImp = y.TDampGroundImp
    x.TDampGroundBare = y.TDampGroundBare
    x.TDampGroundVeg = y.TDampGroundVeg
    x.TDampTree = y.TDampTree
    x.TDampGroundBuild = y.TDampGroundBuild
end

function canyon_fields(::Type{TempDamp})
    return (:TDampGroundImp, :TDampGroundBare, :TDampGroundVeg, :TDampTree)
end

"""
    MRT{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Mean Radiant Temperature Fields.

# Fields
- `Tmrt`: Mean radiant temperature [°C]
- `BoleanInSun`: Point of Tmrt calculation is in sun or in shade
- `SWRdir_Person`: Direct shortwave radiation the person receives [W/m²]
- `SWRdir_in_top`: Direct shortwave radiation the person receives from the top [W/m²]
- `SWRdir_in_bottom`: Direct shortwave radiation the person receives from the bottom [W/m²]
- `SWRdir_in_east`: Direct shortwave radiation the person receives from the east [W/m²]
- `SWRdir_in_south`: Direct shortwave radiation the person receives from the south [W/m²]
- `SWRdir_in_west`: Direct shortwave radiation the person receives from the west [W/m²]
- `SWRdir_in_north`: Direct shortwave radiation the person receives from the north [W/m²]
- `SWRdiff_Person`: Diffuse shortwave radiation the person receives [W/m²]
- `LWR_Person`: Longwave radiation the person receives [W/m²]
"""
Base.@kwdef mutable struct MRT{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    Tmrt::FT
    BoleanInSun::FT
    SWRdir_Person::FT
    SWRdir_in_top::FT
    SWRdir_in_bottom::FT
    SWRdir_in_east::FT
    SWRdir_in_south::FT
    SWRdir_in_west::FT
    SWRdir_in_north::FT
    SWRdiff_Person::FT
    LWR_Person::FT
end

function MRT(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, MRT, Dict{String,Any}())
end

"""
    ThermalComfort{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Universal Thermal Climate Index.

# Fields
- `UTCI`: Universal Thermal Climate Index [°C]
"""
Base.@kwdef mutable struct ThermalComfort{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    UTCI::FT
end

function ThermalComfort(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, ThermalComfort, Dict{String,Any}())
end

"""
    TemperatureVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}

Temperature-related variables for the urban environment.
"""
Base.@kwdef mutable struct TemperatureVariables{FT<:AbstractFloat} <:
                           AbstractModelVariableSet{FT}
    tempvec::TempVec{FT}
    tempdamp::TempDamp{FT}
    mrt::MRT{FT}
    thermalcomfort::ThermalComfort{FT}
end

function TemperatureVariables(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, TemperatureVariables, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{TemperatureVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["tempvec"] = TempVec(FT)
    processed["tempdamp"] = TempDamp(FT)
    processed["mrt"] = MRT(FT)
    processed["thermalcomfort"] = ThermalComfort(FT)
    return processed
end

function ModelComponents.outputs_to_save(
    ::Type{TemperatureVariables}, ::Type{EssentialOutputs}
)
    return (:tempvec, :mrt, :thermalcomfort)
end

function ModelComponents.outputs_to_save(
    ::Type{TemperatureVariables}, ::Type{ExtendedEnergyClimateOutputs}
)
    return (:tempdamp,)
end

function update!(
    x::TemperatureVariables{FT}, results::NamedTuple, ::EBWBCanyonDispatcher
) where {FT<:AbstractFloat}
    _update!(x.tempdamp, results, canyon_fields(TempDamp))
    return nothing
end
