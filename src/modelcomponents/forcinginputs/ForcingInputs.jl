module ForcingInputs

using TethysChlorisCore
using TethysChlorisCore.ModelComponents
using Dates
using NCDatasets

using ..Parameters: LocationProperties
import ...UrbanTethysChloris
using ...UrbanTethysChloris.ModelComponents

abstract type Abstract1PForcingInputs{FT<:AbstractFloat,N} <: AbstractForcingInputs{FT} end
abstract type Abstract1PForcingInputsSet{FT<:AbstractFloat,N} <: AbstractForcingInputSet{FT} end

function Base.getproperty(
    obj::T, field::Symbol
) where {FT<:AbstractFloat,T<:Abstract1PForcingInputs{FT,0}}
    return getfield(obj, field)[]
end

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:Abstract1PForcingInputs{FT,1}}
    scalar_type = typeof(obj).name.wrapper{FT,0}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        fieldvals[field] = fill(getproperty(obj, field)[idx])
    end
    return scalar_type(; fieldvals...)
end

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:Abstract1PForcingInputsSet{FT,1}}
    scalar_type = typeof(obj).name.wrapper{FT,0}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        fieldvals[field] = getindex(getproperty(obj, field), idx)
    end
    return scalar_type(; fieldvals...)
end

include("MeteorologicalInputs.jl")
export MeteorologicalInputs
include("AnthropogenicInputs.jl")
export AnthropogenicInputs
include("SunPositionInputs.jl")
export SunPositionInputs
include("HVACSchedule.jl")
export HVACSchedule
include("ForcingInputSet.jl")

export ForcingInputSet

end
