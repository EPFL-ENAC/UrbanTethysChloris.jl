"""
    SunPositionInputs{FT<:AbstractFloat, N} <: Abstract1PForcingInputs{FT, N}

Solar position and timing parameters.

# Fields
- `t_bef::FT`: Time before solar noon [h]
- `t_aft::FT`: Time after solar noon [h]
- `theta_Z::Vector{FT}`: Solar zenith angle [rad]
- `theta_n::Vector{FT}`: Difference between solar azimuth angle and canyon orientation [rad]
- `zeta_S::Vector{FT}`: Solar azimuth angle [rad]
- `TimeOfMaxSolAlt::Vector{FT}`: Time of maximum solar altitude [h]
"""
Base.@kwdef struct SunPositionInputs{FT<:AbstractFloat,N} <: Abstract1PForcingInputs{FT,N}
    t_bef::FT
    t_aft::FT
    theta_Z::Array{FT,N}
    theta_n::Array{FT,N}
    zeta_S::Array{FT,N}
    TimeOfMaxSolAlt::Array{FT,N}
end

const SCALAR_SUNPOS_FIELDS = [:t_bef, :t_aft]

function TethysChlorisCore.get_required_fields(::Type{SunPositionInputs})
    return [:t_bef, :t_aft]
end

function TethysChlorisCore.get_calculated_fields(::Type{SunPositionInputs})
    return [:theta_Z, :theta_n, :zeta_S, :TimeOfMaxSolAlt]
end

function SunPositionInputs(
    ::Type{FT},
    ::TimeSeries,
    data::NCDataset,
    datetime::Vector{DateTime},
    location::LocationProperties{FT},
) where {FT<:AbstractFloat}
    return initialize(
        FT, SunPositionInputs, data, (FT, dimension_value(TimeSeries())), datetime, location
    )
end

function SunPositionInputs(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return SunPositionInputs{FT,0}(;
        t_bef=FT(data["t_bef"]),
        t_aft=FT(data["t_aft"]),
        theta_Z=fill(FT(data["theta_Z"]), ()),
        theta_n=fill(FT(data["theta_n"]), ()),
        zeta_S=fill(FT(data["zeta_S"]), ()),
        TimeOfMaxSolAlt=fill(FT(data["TimeOfMaxSolAlt"]), ()),
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{SunPositionInputs},
    data::NCDataset,
    params::Tuple,
    datetime::Vector{DateTime},
    location::LocationProperties{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Calculate sun variables for each timestamp
    t_bef = FT(data["t_bef"][])
    t_aft = FT(data["t_aft"][])
    processed["t_bef"] = t_bef
    processed["t_aft"] = t_aft

    results = UrbanTethysChloris.set_sun_variables.(
        datetime, location.DeltaGMT, location.lambda, location.phi, t_bef, t_aft
    )

    h_S = getindex.(results, 1)
    theta_Z = pi/2 .- h_S
    theta_Z[abs.(theta_Z) .≥ π / 2] .= π/2
    processed["theta_Z"] = theta_Z

    zeta_S = getindex.(results, 2)
    processed["zeta_S"] = zeta_S
    processed["theta_n"] = zeta_S .- location.theta_canyon

    T_sunrise = getindex.(results, 3)
    T_sunset = getindex.(results, 4)
    processed["TimeOfMaxSolAlt"] = (T_sunrise .+ T_sunset) ./ 2

    return processed
end

function TethysChlorisCore.validate_fields(::Type{SunPositionInputs}, data::NCDataset) end

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:SunPositionInputs{FT,1}}
    scalar_type = typeof(obj).name.wrapper{FT,0}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        if field in SCALAR_SUNPOS_FIELDS
            fieldvals[field] = getproperty(obj, field)
        else
            fieldvals[field] = fill(getproperty(obj, field)[idx])
        end
    end
    return scalar_type(; fieldvals...)
end
