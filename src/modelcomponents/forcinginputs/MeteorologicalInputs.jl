abstract type AbstractMeteorologicalInputs{FT<:AbstractFloat} <: AbstractForcingInputs{FT} end

Base.@kwdef struct MeteorologicalInputs{FT<:AbstractFloat} <:
                   AbstractMeteorologicalInputs{FT}
    LWR_in::Vector{FT}
    SAB1_in::Vector{FT}
    SAB2_in::Vector{FT}
    SAD1_in::Vector{FT}
    SAD2_in::Vector{FT}
    Tatm::Vector{FT}
    Uatm::Vector{FT}
    Pre::Vector{FT}
    Rain::Vector{FT}
    rel_hum::Vector{FT}
    datetime::Vector{DateTime}
    esat_Tatm::Vector{FT}
    ea::Vector{FT}
    q_atm::Vector{FT}
    qSat_atm::Vector{FT}
    SW_dir::Vector{FT}
    SW_diff::Vector{FT}
    Zatm::FT
    Catm_CO2::FT
    Catm_O2::FT
    SunDSM_MRT::FT
end

function get_scalar_fields(::Type{MeteorologicalInputs})
    return [:Zatm, :Catm_CO2, :Catm_O2, :SunDSM_MRT]
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
        if field ∈ get_scalar_fields(MeteorologicalInputs)
            processed[String(field)] = data[field][]
        else
            processed[String(field)] = Array(data[field])
        end
    end

    processed["Uatm"][processed["Uatm"] .== 0] .= 0.01

    # Calculate ea and q_atm
    processed["esat_Tatm"], processed["ea"] = vapor_pressure(
        processed["Tatm"], processed["rel_hum"]
    )
    processed["q_atm"] = specific_humidity(processed["ea"], processed["Pre"])
    processed["qSat_atm"] = specific_humidity(processed["esat_Tatm"], processed["Pre"])

    processed["SW_dir"] = processed["SAB1_in"] + processed["SAD1_in"]
    processed["SW_diff"] = processed["SAB2_in"] + processed["SAD2_in"]

    update_SW = abs.(cos.(theta_Z)) .< 0.1
    processed["SW_diff"][update_SW] += processed["SW_dir"][update_SW]
    processed["SW_dir"][update_SW] .= 0.0

    return processed
end

function vapor_pressure(Tatm::Vector{<:AbstractFloat}, rel_hum::Vector{<:AbstractFloat})
    esat_Tatm = @. 611.0 * exp(17.27 * (Tatm - 273.16) / (237.3 + Tatm - 273.16))
    ea = esat_Tatm .* rel_hum
    return esat_Tatm, ea
end

function specific_humidity(ea::Vector{<:AbstractFloat}, Pre::Vector{<:AbstractFloat})
    return 0.622 * ea ./ (Pre .- 0.378 .* ea)
end

function TethysChlorisCore.validate_fields(::Type{MeteorologicalInputs}, data::NCDataset) end
