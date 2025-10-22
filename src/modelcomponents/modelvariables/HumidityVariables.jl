abstract type AbstractHumidityVariables{FT<:AbstractFloat,N} <: AbstractModelVariableSet{FT} end

abstract type AbstractHumidityVariablesSubset{FT<:AbstractFloat,N} <:
              AbstractModelVariables{FT} end

abstract type AbstractHumidity{FT<:AbstractFloat,N} <: AbstractHumidityVariablesSubset{FT,N} end
abstract type AbstractResults2m{FT<:AbstractFloat,N} <:
              AbstractHumidityVariablesSubset{FT,N} end

"""
    Humidity{FT<:AbstractFloat, N} <: AbstractHumidity{FT}

Canyon and atmospheric humidity variables.

# Fields
- `CanyonRelative`: Relative humidity at canyon calculation height [-]
- `CanyonSpecific`: Specific humidity at canyon calculation height [kg/kg]
- `CanyonVapourPre`: Vapour pressure at canyon calculation height [Pa]
- `CanyonRelativeSat`: Saturation relative humidity at canyon calculation height [-], is always 1
- `CanyonSpecificSat`: Specific humidity at saturation at canyon calculation height [kg/kg]
- `CanyonVapourPreSat`: Saturation vapor pressure at canyon calculation height [Pa]
- `AtmRelative`: Relative humidity at atmospheric forcing height [-]
- `AtmSpecific`: Specific humidity at atmospheric forcing height [kg/kg]
- `AtmVapourPre`: Vapor pressure at atmospheric forcing height [Pa]
- `AtmRelativeSat`: Saturation relative humidity at atmospheric forcing height [-], is always 1
- `AtmSpecificSat`: Specific humidity at saturation at atmospheric forcing height [kg/kg]
- `AtmVapourPreSat`: Saturation vapour pressure at atmospheric forcing height [Pa]
"""
Base.@kwdef struct Humidity{FT<:AbstractFloat,N} <: AbstractHumidity{FT,N}
    # Canyon and Atmospheric Humidity Fields
    CanyonRelative::Array{FT,N}
    CanyonSpecific::Array{FT,N}
    CanyonVapourPre::Array{FT,N}
    CanyonRelativeSat::Array{FT,N}
    CanyonSpecificSat::Array{FT,N}
    CanyonVapourPreSat::Array{FT,N}
    AtmRelative::Array{FT,N}
    AtmSpecific::Array{FT,N}
    AtmVapourPre::Array{FT,N}
    AtmRelativeSat::Array{FT,N}
    AtmSpecificSat::Array{FT,N}
    AtmVapourPreSat::Array{FT,N}
end

"""
    Results2m{FT<:AbstractFloat, N} <: AbstractResults2m{FT,N}

Temperature and humidity results at 2-meter canyon height.

# Fields
- `T2m`: 2m air temperature [K]
- `q2m`: 2m specific humidity [kg/kg]
- `e_T2m`: 2m vapor pressure [Pa]
- `RH_T2m`: 2m relative humidity [-]
- `qcan`: Canyon specific humidity [kg/kg]
- `e_Tcan`: Canyon vapor pressure [Pa]
- `RH_Tcan`: Canyon relative humidity [-]
"""
Base.@kwdef struct Results2m{FT<:AbstractFloat,N} <: AbstractResults2m{FT,N}
    # 2-meter Height Results
    T2m::Array{FT,N}
    q2m::Array{FT,N}
    e_T2m::Array{FT,N}
    RH_T2m::Array{FT,N}
    qcan::Array{FT,N}
    e_Tcan::Array{FT,N}
    RH_Tcan::Array{FT,N}
end

function Base.getproperty(
    obj::T, field::Symbol
) where {FT<:AbstractFloat,T<:AbstractHumidityVariablesSubset{FT,0}}
    return getfield(obj, field)[]
end

function initialize_humidity(
    ::Type{FT}, N::Int, hours::Int=1, AtmSpecific::FT=0.0
) where {FT<:AbstractFloat}
    return initialize(FT, Humidity, Dict{String,Any}(), (FT, N), hours, AtmSpecific)
end

function initialize_results2m(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, Results2m, Dict{String,Any}(), (FT, N), hours)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    AtmSpecific::FT,
) where {FT<:AbstractFloat,T<:AbstractHumidity}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(T, data, params, hours)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    # Initialize with atmospheric specific humidity if provided
    if params[2] == 1
        processed["CanyonSpecific"][1] = AtmSpecific
    end

    return processed
end

"""
    HumidityVariables{FT<:AbstractFloat, N} <: AbstractHumidityVariables{FT}

    Combines Humidity and Results2m into a single model variables struct.
# Fields
- `Humidity`: Humidity variables at canyon and atmospheric heights
- `Results2m`: Temperature and humidity results at 2-meter canyon height
"""
Base.@kwdef struct HumidityVariables{FT<:AbstractFloat,N} <: AbstractHumidityVariables{FT,N}
    Humidity::Humidity{FT,N}
    Results2m::Results2m{FT,N}
end

function initialize_humidity_variables(
    ::Type{FT}, N::Int, hours::Int=1, AtmSpecific::FT=0.0
) where {FT<:AbstractFloat}
    return initialize(
        FT, HumidityVariables, Dict{String,Any}(), (FT, N), hours, AtmSpecific
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    AtmSpecific::FT,
) where {FT<:AbstractFloat,T<:AbstractHumidityVariables}
    processed = Dict{String,Any}()

    # Preprocess each component separately
    processed["Humidity"] = initialize_humidity(FT, params[2], hours, AtmSpecific)
    processed["Results2m"] = initialize_results2m(FT, params[2], hours)

    return processed
end
