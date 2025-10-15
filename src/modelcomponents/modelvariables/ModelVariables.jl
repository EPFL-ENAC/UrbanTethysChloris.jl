module ModelVariables

using TethysChlorisCore
using TethysChlorisCore.ModelComponents

abstract type AbstractModelVariables{FT<:AbstractFloat} <:
              AbstractIndividualModelComponent{FT} end

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

include("SolverVariables.jl")
include("EnvironmentalConditions.jl")

end
