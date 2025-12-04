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

function Base.show(io::IO, obj::AbstractModelVariables)
    print(io, typeof(obj))
    for field in fieldnames(typeof(obj))
        print(io, "\n", field, ": ", getfield(obj, field))
    end
end

include("BuildingEnergyModelVariables.jl")
export TempVecB, HumidityBuilding, HbuildInt, LEbuildInt, GbuildInt, SWRabsB, LWRabsB
export BEMWasteHeat, BEMEnergyUse, ParACHeat_ts, BuildingEnergyModelVariables
include("EnergyBalanceVariables.jl")
export WBRoof, WBCanyonIndv, WBCanyonTot, EB, SolverVariables, EnergyBalanceVariables
include("EnvironmentalConditions.jl")
export Wind, LAITimeSeries, Resistance, EnvironmentalConditions
include("HeatFluxVariables.jl")
export Hflux, LEflux, Gflux, dStorage, Results2mEnergyFluxes, HeatFluxVariables
include("HumidityVariables.jl")
export Humidity, Results2m, HumidityVariables
include("RadiationFluxVariables.jl")
export AbsorbedRadiationFluxVariablesSubset,
    DefaultRadiationFluxVariablesSubset, AlbedoOutput, RadiationFluxVariables
include("TemperatureVariables.jl")
export TempVec, TempDamp, MRT, ThermalComfort, TemperatureVariables
include("WaterFluxVariables.jl")
export Eflux, Runoff, Runon, Leakage, Interception, dInt_dt, Infiltration, Vwater
export dVwater_dt, Owater, OSwater, Qinlat, ExWater, SoilPotW, CiCO2Leaf, WaterFluxVariables
include("ModelVariableSet.jl")

export Humidity, TempVec, TempVecB
export ModelVariableSet

end
