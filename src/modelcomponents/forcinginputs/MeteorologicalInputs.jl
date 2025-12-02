"""
    MeteorologicalInputs{FT<:AbstractFloat, N}

Meteorological inputs and derived atmospheric properties.

# Fields
- `LWR_in::Vector{FT}`: Atmospheric longwave radiation [W/m² horizontal surface]
- `SAB1_in::Vector{FT}`: Component of direct incoming shortwave radiation [W/m² horizontal surface]
- `SAB2_in::Vector{FT}`: Component of direct incoming shortwave radiation [W/m² horizontal surface]
- `SAD1_in::Vector{FT}`: Component of diffuse incoming shortwave radiation [W/m² horizontal surface]
- `SAD2_in::Vector{FT}`: Component of diffuse incoming shortwave radiation [W/m² horizontal surface]
- `Tatm::Vector{FT}`: Air Temperature at atmospheric reference level [K]
- `Uatm::Vector{FT}`: Wind speed at atmospheric reference level [m/s]
- `Pre::Vector{FT}`: Air pressure [Pa]
- `Rain::Vector{FT}`: Precipitation [mm]
- `rel_hum::Vector{FT}`: Relative humidity [-]
- `datetime::Vector{DateTime}`: Timestamps for the data
- `esat_Tatm::Vector{FT}`: Vapor pressure saturation at Tatm [Pa]
- `ea::Vector{FT}`: Vapor pressure [Pa]
- `q_atm::Vector{FT}`: Specific humidity of air at reference height [-]
- `qSat_atm::Vector{FT}`: Saturation specific humidity [-]
- `SW_dir::Vector{FT}`: Direct incoming shortwave radiation [W/m² horizontal surface]
- `SW_diff::Vector{FT}`: Diffuse incoming shortwave radiation [W/m² horizontal surface]
- `Zatm::FT`: Atmospheric reference height [m]
- `Catm_CO2::FT`: Atmospheric CO2 concentration [ppm]
- `Catm_O2::FT`: Intercellular Partial Pressure Oxygen [ppm]
- `SunDSM_MRT::FT`: Mean radiant temperature from sun [K]
- `cp_atm::Vector{FT}`: Specific heat of air [J/kg K]
- `rho_atm::Vector{FT}`: Dry air density at atmosphere [kg/m³]
"""
Base.@kwdef struct MeteorologicalInputs{FT<:AbstractFloat,N} <:
                   Abstract1PForcingInputs{FT,N}
    LWR_in::Array{FT,N}
    SAB1_in::Array{FT,N}
    SAB2_in::Array{FT,N}
    SAD1_in::Array{FT,N}
    SAD2_in::Array{FT,N}
    Tatm::Array{FT,N}
    Uatm::Array{FT,N}
    Pre::Array{FT,N}
    Rain::Array{FT,N}
    rel_hum::Array{FT,N}
    datetime::Array{DateTime,N}
    esat_Tatm::Array{FT,N}
    ea::Array{FT,N}
    q_atm::Array{FT,N}
    qSat_atm::Array{FT,N}
    SW_dir::Array{FT,N}
    SW_diff::Array{FT,N}
    LWR::Array{FT,N}
    Zatm::FT
    Catm_CO2::FT
    Catm_O2::FT
    SunDSM_MRT::FT
    cp_atm::Array{FT,N}
    rho_atm::Array{FT,N}
    AtmRelative::Array{FT,N}
    AtmSpecific::Array{FT,N}
    AtmVapourPre::Array{FT,N}
    AtmRelativeSat::Array{FT,N}
    AtmSpecificSat::Array{FT,N}
    AtmVapourPreSat::Array{FT,N}
end

const SCALAR_METEO_FIELDS = [:Zatm, :Catm_CO2, :Catm_O2, :SunDSM_MRT]

function get_scalar_fields(::Type{MeteorologicalInputs})
    return [:Zatm, :Catm_CO2, :Catm_O2, :SunDSM_MRT]
end

function TethysChlorisCore.get_calculated_fields(::Type{MeteorologicalInputs})
    return [
        :esat_Tatm,
        :ea,
        :q_atm,
        :qSat_atm,
        :SW_dir,
        :SW_diff,
        :cp_atm,
        :rho_atm,
        :AtmRelative,
        :AtmSpecific,
        :AtmVapourPre,
        :AtmRelativeSat,
        :AtmSpecificSat,
        :AtmVapourPreSat,
        :LWR,
    ]
end

function MeteorologicalInputs(
    ::Type{FT}, ::TimeSeries, data::NCDataset, theta_Z::Vector{FT}
) where {FT<:AbstractFloat}
    return initialize(
        FT, MeteorologicalInputs, data, (FT, dimension_value(TimeSeries())), theta_Z
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{MeteorologicalInputs},
    data::NCDataset,
    params::Tuple,
    theta_Z::Vector{FT},
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

    processed["cp_atm"] =
        1005.0 .+ (((processed["Tatm"] .- 273.15) .+ 23.15) .^ 2) .+ 3364.0
    processed["rho_atm"] =
        (processed["Pre"] ./ (287.04 .* processed["Tatm"])) .*
        (1.0 .- (processed["ea"] ./ processed["Pre"]) .* (1.0 - 0.622))

    processed["AtmRelative"] = processed["rel_hum"]
    processed["AtmSpecific"] = processed["q_atm"]
    processed["AtmVapourPre"] = processed["ea"]
    processed["AtmRelativeSat"] = ones(size(processed["ea"]))
    processed["AtmSpecificSat"] = processed["qSat_atm"]
    processed["AtmVapourPreSat"] = processed["esat_Tatm"]
    processed["LWR"] = processed["LWR_in"]

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

function Base.getindex(
    obj::T, idx::Int
) where {FT<:AbstractFloat,T<:MeteorologicalInputs{FT,1}}
    scalar_type = typeof(obj).name.wrapper{FT,0}
    fieldvals = Dict{Symbol,Any}()
    for field in fieldnames(typeof(obj))
        if field in SCALAR_METEO_FIELDS
            fieldvals[field] = getproperty(obj, field)
        else
            fieldvals[field] = fill(getproperty(obj, field)[idx])
        end
    end
    return scalar_type(; fieldvals...)
end
