"""
    Humidity{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct Humidity{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
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

function initialize_humidity(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Humidity, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_humidity(
    ::Type{FT}, ::TimeSeries, hours::Int, AtmSpecific::FT
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        Humidity,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
        AtmSpecific,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{Humidity},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    AtmSpecific::FT,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(Humidity, data, params, hours)

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
Base.@kwdef mutable struct Results2m{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    # 2-meter Height Results
    T2m::Array{FT,N}
    q2m::Array{FT,N}
    e_T2m::Array{FT,N}
    RH_T2m::Array{FT,N}
    qcan::Array{FT,N}
    e_Tcan::Array{FT,N}
    RH_Tcan::Array{FT,N}
end

function initialize_results2m(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Results2m, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_results2m(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT, Results2m, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    HumidityVariables{FT<:AbstractFloat, N} <: Abstract1PModelVariablesSet{FT, N}

    Combines Humidity and Results2m into a single model variables struct.
# Fields
- `Humidity`: Humidity variables at canyon and atmospheric heights
- `Results2m`: Temperature and humidity results at 2-meter canyon height
"""
Base.@kwdef struct HumidityVariables{FT<:AbstractFloat,N} <:
                   Abstract1PModelVariablesSet{FT,N}
    Humidity::Humidity{FT,N}
    Results2m::Results2m{FT,N}
end

function initialize_humidity_variables(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, HumidityVariables, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{HumidityVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Preprocess each component separately
    processed["Humidity"] = initialize_humidity(FT, dimensionality_type(params[2]))
    processed["Results2m"] = initialize_results2m(FT, dimensionality_type(params[2]))

    return processed
end

function initialize_humidity_variables(
    ::Type{FT}, ::TimeSeries, hours::Int, AtmSpecific::FT
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        HumidityVariables,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
        AtmSpecific,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{HumidityVariables},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    AtmSpecific::FT,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Preprocess each component separately
    processed["Humidity"] = initialize_humidity(
        FT, dimensionality_type(params[2]), hours, AtmSpecific
    )
    processed["Results2m"] = initialize_results2m(FT, dimensionality_type(params[2]), hours)

    return processed
end
