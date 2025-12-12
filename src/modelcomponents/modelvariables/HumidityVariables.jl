"""
    Humidity{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct Humidity{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    CanyonRelative::FT
    CanyonSpecific::FT
    CanyonVapourPre::FT
    CanyonRelativeSat::FT
    CanyonSpecificSat::FT
    CanyonVapourPreSat::FT
    AtmRelative::FT
    AtmSpecific::FT
    AtmVapourPre::FT
    AtmRelativeSat::FT
    AtmSpecificSat::FT
    AtmVapourPreSat::FT
    q2m::FT
end

function Humidity(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Humidity, Dict{String,Any}())
end

function Humidity(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return Humidity{FT}(;
        CanyonRelative=data["CanyonRelative"],
        CanyonSpecific=data["CanyonSpecific"],
        CanyonVapourPre=data["CanyonVapourPre"],
        CanyonRelativeSat=data["CanyonRelativeSat"],
        CanyonSpecificSat=data["CanyonSpecificSat"],
        CanyonVapourPreSat=data["CanyonVapourPreSat"],
        AtmRelative=data["AtmRelative"],
        AtmSpecific=data["AtmSpecific"],
        AtmVapourPre=data["AtmVapourPre"],
        AtmRelativeSat=data["AtmRelativeSat"],
        AtmSpecificSat=data["AtmSpecificSat"],
        AtmVapourPreSat=data["AtmVapourPreSat"],
        q2m=data["q2m"],
    )
end

function update!(dest::Humidity{FT}, src::Humidity{FT}) where {FT<:AbstractFloat}
    dest.CanyonRelative = src.CanyonRelative
    dest.CanyonSpecific = src.CanyonSpecific
    dest.CanyonVapourPre = src.CanyonVapourPre
    dest.CanyonRelativeSat = src.CanyonRelativeSat
    dest.CanyonSpecificSat = src.CanyonSpecificSat
    dest.CanyonVapourPreSat = src.CanyonVapourPreSat
    dest.AtmRelative = src.AtmRelative
    dest.AtmSpecific = src.AtmSpecific
    dest.AtmVapourPre = src.AtmVapourPre
    dest.AtmRelativeSat = src.AtmRelativeSat
    dest.AtmSpecificSat = src.AtmSpecificSat
    dest.AtmVapourPreSat = src.AtmVapourPreSat
    dest.q2m = src.q2m
end

"""
    Results2m{FT<:AbstractFloat} <: AbstractResults2m{FT}

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
Base.@kwdef mutable struct Results2m{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    T2m::FT
    q2m::FT
    e_T2m::FT
    RH_T2m::FT
    qcan::FT
    e_Tcan::FT
    RH_Tcan::FT
end

function Results2m(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Results2m, Dict{String,Any}())
end

"""
    HumidityVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}

    Combines Humidity and Results2m into a single model variables struct.
# Fields
- `Humidity`: Humidity variables at canyon and atmospheric heights
- `Results2m`: Temperature and humidity results at 2-meter canyon height
"""
Base.@kwdef struct HumidityVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}
    Humidity::Humidity{FT}
    Results2m::Results2m{FT}
end

function HumidityVariables(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, HumidityVariables, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{HumidityVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Preprocess each component separately
    processed["Humidity"] = Humidity(FT)
    processed["Results2m"] = Results2m(FT)

    return processed
end
