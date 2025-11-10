module ModelVariables

using TethysChlorisCore
using TethysChlorisCore.ModelComponents
using ..Parameters:
    SoilParameters,
    VegetationParameters,
    HeightDependentVegetationParameters,
    VegetatedSoilParameters
using ...Soil: soil_parameters_total

abstract type AbstractModelVariables{FT<:AbstractFloat} <:
              AbstractIndividualModelComponent{FT} end
abstract type Abstract1PModelVariables{FT<:AbstractFloat,N} <: AbstractModelVariables{FT} end
abstract type Abstract2PModelVariables{FT<:AbstractFloat,N,M} <: AbstractModelVariables{FT} end

abstract type AbstractModelVariableSet{FT<:AbstractFloat} <: AbstractModelComponentSet{FT} end
abstract type Abstract1PModelVariablesSet{FT<:AbstractFloat,N} <:
              AbstractModelVariableSet{FT} end
abstract type Abstract2PModelVariablesSet{FT<:AbstractFloat,N,Np1} <:
              AbstractModelVariableSet{FT} end

abstract type ModelDimension end
struct TimeSlice <: ModelDimension end
struct TimeSeries <: ModelDimension end
dimension_value(::TimeSlice) = 0
dimension_value(::TimeSeries) = 1
dimensionality_type(dim_value::Int) = dim_value == 0 ? TimeSlice() : TimeSeries()

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{T}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat,T<:AbstractModelVariables}
    processed = Dict{String,Any}()

    vars = String.(fieldnames(T))

    for var in vars
        processed[var] = zeros(FT, ())
    end

    return processed
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{T}, data::Dict{String,Any}, params::Tuple, args...
) where {FT<:AbstractFloat,T<:AbstractModelVariables}
    processed = Dict{String,Any}()

    dimensions = get_dimensions(T, data, params, args...)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    return processed
end

function TethysChlorisCore.get_required_fields(
    ::Type{T}
) where {T<:Union{AbstractModelVariables,AbstractModelVariableSet}}
    return []
end

"""
    get_dimensions(
        ::Type{T}, data::Dict{String,Any}, params::Tuple{DataType, Signed},
    ) where {T<:AbstractModelVariables}

Get the dimensions of each field in the model variables struct `T` based on the parameters.
"""
function get_dimensions(
    ::Type{T}, data::Dict{String,Any}, params::Tuple{DataType,Signed}
) where {T<:AbstractModelVariables}
    dimensions = Dict{String,Tuple}()

    for field in fieldnames(T)
        dimensions[String(field)] = ()
    end

    return dimensions
end

"""
    get_dimensions(
        ::Type{T}, data::Dict{String,Any}, params::Tuple{DataType, Signed}, hours::Int
    ) where {T<:AbstractModelVariables}

Get the dimensions of each field in the model variables struct `T` based on the parameters and hours.
"""
function get_dimensions(
    ::Type{T}, data::Dict{String,Any}, params::Tuple{DataType,Signed}, hours::Int
) where {T<:AbstractModelVariables}
    dimensions = Dict{String,Tuple}()

    for field in fieldnames(T)
        dimensions[String(field)] = (hours,)
    end

    return dimensions
end

# getproperty specializations for model variables
function Base.getproperty(
    obj::T, field::Symbol
) where {FT<:AbstractFloat,T<:Abstract1PModelVariables{FT,0}}
    return getfield(obj, field)[]
end

function Base.getproperty(
    obj::T, field::Symbol
) where {FT<:AbstractFloat,T<:Abstract2PModelVariables{FT,0,1}}
    if field in get_vector_fields(obj)
        return getfield(obj, field)[]
    else
        return getfield(obj, field)
    end
end

# setproperty! specializations for model variables
function Base.setproperty!(
    obj::T, field::Symbol, value::FT
) where {FT<:AbstractFloat,T<:Abstract1PModelVariables{FT,0}}
    setfield!(obj, field, fill(value))
end

# getindex specializations for model variables
function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:Abstract1PModelVariables{FT,1}}
    scalar_type = typeof(obj).name.wrapper{FT,0}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        fieldvals[field] = fill(getproperty(obj, field)[idx])
    end
    return scalar_type(; fieldvals...)
end

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:Abstract1PModelVariables{FT,2}}
    scalar_type = typeof(obj).name.wrapper{FT,1}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        fieldvals[field] = getproperty(obj, field)[idx, :]
    end
    return scalar_type(; fieldvals...)
end

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:Abstract2PModelVariables{FT,1,2}}
    scalar_type = typeof(obj).name.wrapper{FT,0,1}
    fieldvals = Dict{Symbol,Any}()

    vector_fields = get_vector_fields(obj)

    for field in fieldnames(typeof(obj))
        if field in vector_fields
            fieldvals[field] = fill(getproperty(obj, field)[idx])
        else
            fieldvals[field] = getproperty(obj, field)[idx, :]
        end
    end
    return scalar_type(; fieldvals...)
end

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:Abstract1PModelVariablesSet{FT,1}}
    scalar_type = typeof(obj).name.wrapper{FT,0}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        fieldvals[field] = getindex(getproperty(obj, field), idx)
    end
    return scalar_type(; fieldvals...)
end

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:Abstract2PModelVariablesSet{FT,1,2}}
    scalar_type = typeof(obj).name.wrapper{FT,0,1}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        fieldvals[field] = getindex(getproperty(obj, field), idx)
    end
    return scalar_type(; fieldvals...)
end

# setindex! specializations for model variables
function Base.setindex!(
    obj::T, value::S, idx::Int
) where {
    FT<:AbstractFloat,T<:Abstract1PModelVariables{FT,1},S<:Abstract1PModelVariables{FT,0}
}
    for field in fieldnames(typeof(obj))
        getproperty(obj, field)[idx] = getproperty(value, field)
    end
    return obj
end

function Base.setindex!(
    obj::T, value::S, idx::Int
) where {
    FT<:AbstractFloat,T<:Abstract1PModelVariables{FT,2},S<:Abstract1PModelVariables{FT,1}
}
    for field in fieldnames(typeof(obj))
        getproperty(obj, field)[idx, :] = getproperty(value, field)
    end
    return obj
end

function Base.setindex!(
    obj::T, value::S, idx::Int
) where {
    FT<:AbstractFloat,
    T<:Abstract2PModelVariables{FT,1,2},
    S<:Abstract2PModelVariables{FT,0,1},
}
    vector_fields = get_vector_fields(obj)

    for field in fieldnames(typeof(obj))
        if field in vector_fields
            getproperty(obj, field)[idx] = getproperty(value, field)
        else
            getproperty(obj, field)[idx, :] = getproperty(value, field)
        end
    end

    return obj
end

function Base.setindex!(
    obj::T, value::S, idx::Int
) where {
    FT<:AbstractFloat,
    T<:Abstract1PModelVariablesSet{FT,1},
    S<:Abstract1PModelVariablesSet{FT,0},
}
    for field in fieldnames(typeof(obj))
        getproperty(obj, field)[idx] = getproperty(value, field)
    end

    return obj
end

function Base.setindex!(
    obj::T, value::S, idx::Int
) where {
    FT<:AbstractFloat,
    T<:Abstract2PModelVariablesSet{FT,1,2},
    S<:Abstract2PModelVariablesSet{FT,0,1},
}
    for field in fieldnames(typeof(obj))
        getproperty(obj, field)[idx] = getproperty(value, field)
    end
    return obj
end

include("BuildingEnergyModelVariables.jl")
include("EnergyBalanceVariables.jl")
include("EnvironmentalConditions.jl")
include("HeatFluxVariables.jl")
include("HumidityVariables.jl")
include("RadiationFluxVariables.jl")
include("TemperatureVariables.jl")
include("WaterFluxVariables.jl")
include("ModelVariableSet.jl")

end
