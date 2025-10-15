module ModelVariables

using TethysChlorisCore
using TethysChlorisCore.ModelComponents

abstract type AbstractModelVariables{FT<:AbstractFloat} <:
              AbstractIndividualModelComponent{FT} end

abstract type AbstractModelVariableSet{FT<:AbstractFloat} <: AbstractModelComponentSet{FT} end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{T}, data::Dict{String,Any}, params::Tuple, hours::Int
) where {FT<:AbstractFloat,T<:AbstractModelVariables}
    processed = Dict{String,Any}()

    dimensions = get_dimensions(T, data, params, hours)

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
        ::Type{T}, data::Dict{String,Any}, params::Tuple{DataType, Signed}, hours::Int
    ) where {T<:AbstractModelVariables}

Get the dimensions of each field in the model variables struct `T` based on the parameters and hours.
"""
function get_dimensions(
    ::Type{T}, data::Dict{String,Any}, params::Tuple{DataType,Signed}, hours::Int
) where {T<:AbstractModelVariables}
    if params[2] ∉ [0, 1]
        throw(ArgumentError("Only N=0 and N=1 are currently supported"))
    end

    # Create a dictionary with all field names and their dimensions
    dimensions = Dict{String,Tuple}()

    # Get all field names from the struct
    field_names = fieldnames(T)

    if params[2] == 0
        # For scalar case (N=0), all fields are scalars
        for field in field_names
            dimensions[String(field)] = ()
        end
    else
        # For time series case (N=1), all fields have hours dimension
        for field in field_names
            dimensions[String(field)] = (hours,)
        end
    end

    return dimensions
end

include("BuildingEnergyModelVariables.jl")
include("EnvironmentalConditions.jl")
include("HeatFluxVariables.jl")
include("HumidityVariables.jl")
include("SolverVariables.jl")
include("TemperatureVariables.jl")

end
