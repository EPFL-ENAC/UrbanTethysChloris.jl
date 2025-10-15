"""
    SolverVariables{FT<:AbstractFloat, N, Np} <: AbstractModelVariables{FT}

Optical properties for indoor building surfaces.

# Fields
- `Success` Boolean indicating convergence of solution of energy balance
- `ValuesEB`: Energy balance closure for the different equations [W/m²]
- `Tsolver`: Temperatures and humidity of different canyon faces and air [K], [kg/kg]
- `YfunctionOutput`: Solver function outputs
"""
Base.@kwdef struct SolverVariables{FT<:AbstractFloat,N,Np} <: AbstractModelVariables{FT}
    Success::Array{Bool,N}
    ValuesEB::Array{FT,Np}
    Tsolver::Array{FT,Np}
    YfunctionOutput::Array{FT,Np}
end

function TethysChlorisCore.get_required_fields(::Type{SolverVariables})
    return []
end

function get_vector_fields(obj::SolverVariables{FT,0,1}) where {FT<:AbstractFloat}
    return [:Success]
end

function Base.getproperty(
    obj::SolverVariables{FT,0,1}, field::Symbol
) where {FT<:AbstractFloat}
    if field in get_vector_fields(obj)
        return getfield(obj, field)[]
    else
        return getfield(obj, field)
    end
end

function initialize_solver_variables(
    ::Type{FT}, N::Int, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, SolverVariables, Dict{String,Any}(), (FT, N, N+1), hours)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{T}, data::Dict{String,Any}, params::Tuple, hours::Int
) where {FT<:AbstractFloat,T<:AbstractModelVariables}
    processed = Dict{String,Any}()

    dimensions = get_dimensions(T, data, params, hours)

    for (var, dims) in dimensions
        if var != "Success"
            processed[var] = zeros(FT, dims)
        else
            processed[var] = fill(false, dims)
        end
    end

    return processed
end

function get_dimensions(
    ::Type{SolverVariables}, data::Dict{String,Any}, params::Tuple, hours::Int
)
    if params[2] ∉ [0, 1]
        throw(ArgumentError("Only N=0 and N=1 are currently supported"))
    end

    if params[2] == 0
        return Dict(
            "Success" => (),
            "ValuesEB" => (22,),
            "Tsolver" => (22,),
            "YfunctionOutput" => (22,),
        )
    else
        return Dict(
            "Success" => (hours,),
            "ValuesEB" => (hours, 22),
            "Tsolver" => (hours, 22),
            "YfunctionOutput" => (hours, 22),
        )
    end
end
