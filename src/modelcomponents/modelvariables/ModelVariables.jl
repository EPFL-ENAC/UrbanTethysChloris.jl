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

abstract type AbstractModelVariableSet{FT<:AbstractFloat} <: AbstractModelComponentSet{FT} end

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

Get the dimensions of each field in the model variables struct `T` based on the parameters and hours.
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
