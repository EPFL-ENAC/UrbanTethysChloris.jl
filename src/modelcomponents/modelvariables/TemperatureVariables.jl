"""
    TempVec{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Temperature Vector Fields.

# Fields
- `TRoofImp`: Temperature roof impervious area [K]
- `TRoofVeg`: Temperature roof vegetated area [K]
- `TRoofIntImp`: Interior temperature roof impervious area [K]
- `TRoofIntVeg`: Interior temperature roof vegetated area [K]
- `TGroundImp`: Temperature ground impervious area [K]
- `TGroundBare`: Temperature ground bare area [K]
- `TGroundVeg`: Temperature ground vegetated area [K]
- `TTree`: Temperature tree canopy [K]
- `TWallSun`: Temperature sunlit wall area [K]
- `TWallShade`: Temperature shaded wall area [K]
- `TWallIntSun`: Interior temperature sunlit wall [K]
- `TWallIntShade`: Interior temperature shaded wall [K]
- `TCanyon`: Temperature canyon [K]
- `Tatm`: Temperature atmosphere (measured) [K]
"""
Base.@kwdef mutable struct TempVec{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    TRoofImp::Array{FT,N}
    TRoofVeg::Array{FT,N}
    TRoofIntImp::Array{FT,N}
    TRoofIntVeg::Array{FT,N}
    TGroundImp::Array{FT,N}
    TGroundBare::Array{FT,N}
    TGroundVeg::Array{FT,N}
    TTree::Array{FT,N}
    TWallSun::Array{FT,N}
    TWallShade::Array{FT,N}
    TWallIntSun::Array{FT,N}
    TWallIntShade::Array{FT,N}
    TCanyon::Array{FT,N}
    Tatm::Array{FT,N}
end

function initialize_tempvec(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, TempVec, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_tempvec(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, TempVec, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    TempDamp{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Temperature Dampening Fields.

# Fields
- `TDampGroundImp`: Dampening temperature ground impervious area [K]
- `TDampGroundBare`: Dampening temperature ground bare area [K]
- `TDampGroundVeg`: Dampening temperature ground vegetated area [K]
- `TDampTree`: Dampening temperature tree canopy [K]
- `TDampGroundBuild`: Dampening temperature of ground in building interior [K]
"""
Base.@kwdef mutable struct TempDamp{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    TDampGroundImp::Array{FT,N}
    TDampGroundBare::Array{FT,N}
    TDampGroundVeg::Array{FT,N}
    TDampTree::Array{FT,N}
    TDampGroundBuild::Array{FT,N}
end

function initialize_tempdamp(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, TempDamp, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_tempdamp(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, TempDamp, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    MRT{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Mean Radiant Temperature Fields.

# Fields
- `Tmrt`: Mean radiant temperature [°C]
- `BoleanInSun`: Point of Tmrt calculation is in sun or in shade
- `SWRdir_Person`: Direct shortwave radiation the person receives [W/m²]
- `SWRdir_in_top`: Direct shortwave radiation the person receives from the top [W/m²]
- `SWRdir_in_bottom`: Direct shortwave radiation the person receives from the bottom [W/m²]
- `SWRdir_in_east`: Direct shortwave radiation the person receives from the east [W/m²]
- `SWRdir_in_south`: Direct shortwave radiation the person receives from the south [W/m²]
- `SWRdir_in_west`: Direct shortwave radiation the person receives from the west [W/m²]
- `SWRdir_in_north`: Direct shortwave radiation the person receives from the north [W/m²]
- `SWRdiff_Person`: Diffuse shortwave radiation the person receives [W/m²]
- `LWR_Person`: Longwave radiation the person receives [W/m²]
"""
Base.@kwdef mutable struct MRT{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    Tmrt::Array{FT,N}
    BoleanInSun::Array{FT,N}
    SWRdir_Person::Array{FT,N}
    SWRdir_in_top::Array{FT,N}
    SWRdir_in_bottom::Array{FT,N}
    SWRdir_in_east::Array{FT,N}
    SWRdir_in_south::Array{FT,N}
    SWRdir_in_west::Array{FT,N}
    SWRdir_in_north::Array{FT,N}
    SWRdiff_Person::Array{FT,N}
    LWR_Person::Array{FT,N}
end

function initialize_mrt(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, MRT, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_mrt(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, MRT, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    ThermalComfort{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Universal Thermal Climate Index.

# Fields
- `UTCI`: Universal Thermal Climate Index [°C]
"""
Base.@kwdef mutable struct ThermalComfort{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    UTCI::Array{FT,N}
end

function initialize_thermal_comfort(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, ThermalComfort, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_thermal_comfort(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT, ThermalComfort, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    TemperatureVariables{FT<:AbstractFloat} <: Abstract1PModelVariablesSet{FT, N}

Temperature-related variables for the urban environment.
"""
Base.@kwdef mutable struct TemperatureVariables{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariablesSet{FT,N}
    tempvec::TempVec{FT,N}
    tempdamp::TempDamp{FT,N}
    mrt::MRT{FT,N}
    thermalcomfort::ThermalComfort{FT,N}
end

function initialize_temperature_variables(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, TemperatureVariables, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_temperature_variables(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        TemperatureVariables,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{TemperatureVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["tempvec"] = initialize_tempvec(FT, dimensionality_type(params[2]))
    processed["tempdamp"] = initialize_tempdamp(FT, dimensionality_type(params[2]))
    processed["mrt"] = initialize_mrt(FT, dimensionality_type(params[2]))
    processed["thermalcomfort"] = initialize_thermal_comfort(
        FT, dimensionality_type(params[2])
    )
    return processed
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{TemperatureVariables},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["tempvec"] = initialize_tempvec(FT, dimensionality_type(params[2]), hours)
    processed["tempdamp"] = initialize_tempdamp(FT, dimensionality_type(params[2]), hours)
    processed["mrt"] = initialize_mrt(FT, dimensionality_type(params[2]), hours)
    processed["thermalcomfort"] = initialize_thermal_comfort(
        FT, dimensionality_type(params[2]), hours
    )
    return processed
end
