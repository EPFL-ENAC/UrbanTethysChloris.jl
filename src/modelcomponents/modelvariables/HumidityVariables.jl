abstract type AbstractHumidityVariables{FT<:AbstractFloat} <: AbstractModelVariables{FT} end

"""
    HumidityVariables{FT<:AbstractFloat, N} <: AbstractModelVariables{FT}

Humidity-related variables for the urban environment.

# Fields
## Canyon and Atmospheric Humidity Fields
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

## 2-meter Height Results
- `T2m`: 2m air temperature [K]
- `q2m`: 2m specific humidity [kg/kg]
- `e_T2m`: 2m vapor pressure [Pa]
- `RH_T2m`: 2m relative humidity [-]
- `qcan`: Canyon specific humidity [kg/kg]
- `e_Tcan`: Canyon vapor pressure [Pa]
- `RH_Tcan`: Canyon relative humidity [-]
"""
Base.@kwdef struct HumidityVariables{FT<:AbstractFloat,N} <: AbstractHumidityVariables{FT}
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
    obj::HumidityVariables{FT,0}, field::Symbol
) where {FT<:AbstractFloat}
    return getfield(obj, field)[]
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
) where {FT<:AbstractFloat,T<:AbstractModelVariables}
    processed = Dict{String,Any}()

    dimensions = get_dimensions(T, data, params, hours)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    if params[2] == 1
        processed["CanyonSpecific"][1] = AtmSpecific
    end

    return processed
end
