module ModelVariables

using TethysChlorisCore
using TethysChlorisCore.ModelComponents
using ..Parameters:
    SoilParameters,
    VegetationParameters,
    HeightDependentVegetationParameters,
    VegetatedSoilParameters
using ...UrbanTethysChloris.ModelComponents
using StaticArrays

abstract type AbstractModelVariables{FT<:AbstractFloat} <:
              AbstractIndividualModelComponent{FT} end
abstract type AbstractModelVariableSet{FT<:AbstractFloat} <: AbstractModelComponentSet{FT} end
abstract type AbstractLayeredSoilVariables{FT<:AbstractFloat} <: AbstractModelVariables{FT} end

function TethysChlorisCore.get_required_fields(
    ::Type{T}
) where {T<:Union{AbstractModelVariables,AbstractModelVariableSet}}
    return []
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{T}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat,T<:AbstractModelVariables}
    processed = Dict{String,Any}()

    vars = String.(fieldnames(T))

    for var in vars
        processed[var] = zero(FT)
    end

    return processed
end

function allocate_fields!(
    x::Dict{String,Any},
    ::Type{FT},
    ::Type{T},
    fields::Vector{String},
    vector_length::Signed,
) where {FT<:AbstractFloat,T<:AbstractLayeredSoilVariables}
    for field in fields
        x[field] = zeros(FT, vector_length)
    end
end

# Default initialization does nothing
function initialize_fields!(
    x::Dict{String,Any},
    ::Type{T},
    fields::Vector{String},
    soil::VegetatedSoilParameters{FT},
) where {T<:AbstractLayeredSoilVariables,FT<:AbstractFloat} end

function add_fields!(
    x::Dict{String,Any},
    ::Type{FT},
    ::Type{T},
    fields::Vector{String},
    soil::VegetatedSoilParameters{FT},
) where {FT<:AbstractFloat,T<:AbstractLayeredSoilVariables}
    allocate_fields!(x, FT, T, fields, soil.ms)

    initialize_fields!(x, T, fields, soil)
end

function roof_fields(::Type{T}) where {T<:AbstractLayeredSoilVariables}
    return String[]
end

function ground_fields(::Type{T}) where {T<:AbstractLayeredSoilVariables}
    return String[]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    soil::SoilParameters{FT},
    args...,
) where {FT<:AbstractFloat,T<:AbstractLayeredSoilVariables}
    processed = Dict{String,Any}()

    add_fields!(processed, FT, T, ground_fields(T), soil.ground, args...)

    add_fields!(processed, FT, T, roof_fields(T), soil.roof, args...)

    return processed
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

export Humidity, TempVec, TempVecB
export ModelVariableSet

end
