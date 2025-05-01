abstract type AbstractMeteorologicalInputs{FT<:AbstractFloat} <: AbstractForcingInputs{FT} end

Base.@kwdef struct MeteorologicalInputs{FT<:AbstractFloat} <:
                   AbstractMeteorologicalInputs{FT}
    LWR_in::Vector{FT}
    SAB1_in::Vector{FT}
    SAB2_in::Vector{FT}
    SAD1_in::Vector{FT}
    SAD2_in::Vector{FT}
    T_atm::Vector{FT}
    windspeed_u::Vector{FT}
    pressure_atm::Vector{FT}
    rain::Vector{FT}
    rel_humidity::Vector{FT}
    datetime::Vector{DateTime}
    esat_Tatm::Vector{FT}
    ea::Vector{FT}
    q_atm::Vector{FT}
    qSat_atm::Vector{FT}
    SW_dir::Vector{FT}
    SW_diff::Vector{FT}
end

function TethysChlorisCore.get_calculated_fields(::Type{MeteorologicalInputs})
    return [:esat_Tatm, :ea, :q_atm, :qSat_atm, :SW_dir, :SW_diff]
end

function initialize_meteorological_inputs(
    ::Type{FT}, data::NCDataset, theta_Z::Vector{FT}
) where {FT<:AbstractFloat}
    return initialize(FT, MeteorologicalInputs, data, theta_Z)
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{MeteorologicalInputs}, data::NCDataset, theta_Z::Vector{FT}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    fields = [
        TethysChlorisCore.get_required_fields(MeteorologicalInputs);
        filter(
            name -> haskey(data, name),
            TethysChlorisCore.get_optional_fields(MeteorologicalInputs),
        )
    ]

    for field in fields
        processed[String(field)] = data[field][:]
    end

    # Calculate ea and q_atm
    processed["esat_Tatm"], processed["ea"] = vapor_pressure(
        processed["T_atm"], processed["rel_humidity"]
    )
    processed["q_atm"] = specific_humidity(processed["ea"], processed["pressure_atm"])
    processed["qSat_atm"] = specific_humidity(
        processed["esat_Tatm"], processed["pressure_atm"]
    )

    processed["SW_dir"] = processed["SAB1_in"] + processed["SAD1_in"]
    processed["SW_diff"] = processed["SAB2_in"] + processed["SAD2_in"]

    update_SW = abs.(cos.(theta_Z)) .< 0.1
    processed["SW_diff"][update_SW] += processed["SW_dir"][update_SW]
    processed["SW_dir"][update_SW] .= 0.0

    return processed
end

function vapor_pressure(
    Tatm::Vector{<:AbstractFloat}, rel_humidity::Vector{<:AbstractFloat}
)
    esat_Tatm = @. 611.0 * exp(17.27 * (Tatm - 273.16) / (237.3 + Tatm - 273.16))
    ea = esat_Tatm .* rel_humidity
    return esat_Tatm, ea
end

function specific_humidity(
    ea::Vector{<:AbstractFloat}, pressure_atm::Vector{<:AbstractFloat}
)
    return 0.622 * ea ./ (pressure_atm .- 0.378 .* ea)
end

function TethysChlorisCore.validate_fields(::Type{MeteorologicalInputs}, data::NCDataset) end
