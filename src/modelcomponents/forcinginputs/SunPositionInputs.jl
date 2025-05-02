abstract type AbstractSunPositionInputs{FT<:AbstractFloat} <: AbstractForcingInputs{FT} end

Base.@kwdef struct SunPositionInputs{FT<:AbstractFloat} <: AbstractSunPositionInputs{FT}
    t_bef::FT
    t_aft::FT
    theta_Z::Vector{FT}
    theta_n::Vector{FT}
    zeta_S::Vector{FT}
    TimeOfMaxSolAlt::Vector{FT}
end

function TethysChlorisCore.get_required_fields(::Type{SunPositionInputs})
    return [:t_bef, :t_aft]
end

function TethysChlorisCore.get_calculated_fields(::Type{SunPositionInputs})
    return [:theta_Z, :theta_n, :zeta_S, :TimeOfMaxSolAlt]
end

function initialize_sunposition_inputs(
    ::Type{FT},
    data::NCDataset,
    datetime::Vector{DateTime},
    location::LocationProperties{FT},
) where {FT<:AbstractFloat}
    return initialize(FT, SunPositionInputs, data, datetime, location)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{SunPositionInputs},
    data::NCDataset,
    datetime::Vector{DateTime},
    location::LocationProperties{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Calculate sun variables for each timestamp
    t_bef = data["t_bef"][]
    t_aft = data["t_aft"][]
    processed["t_bef"] = t_bef
    processed["t_aft"] = t_aft

    results = UrbanTethysChloris.set_sun_variables.(
        datetime, location.DeltaGMT, location.lambda, location.phi, t_bef, t_aft
    )

    h_S = getindex.(results, 1)
    theta_Z = pi/2 .- h_S
    theta_Z[abs.(theta_Z) .≥ -π / 2] .= π/2
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
