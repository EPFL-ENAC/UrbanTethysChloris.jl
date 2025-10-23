"""
    ModelVariableSet{FT<:AbstractFloat,N,Np} <: AbstractModelVariableSet{FT}

Forcing inputs for the Urban Tethys-Chloris model.

# Fields
- `buildingenergymodel`: Building energy model variables
- `energybalance`: Energy balance variables
- `environmentalconditions`: Environmental conditions variables
- `heatflux`: Heat flux variables
- `humidity`: Humidity variables
- `radiationflux`: Radiation flux variables
- `temperature`: Temperature variables
- `waterflux`: Water flux variables
"""
Base.@kwdef struct ModelVariableSet{FT<:AbstractFloat,N,Np} <: AbstractModelVariableSet{FT}
    buildingenergymodel::BuildingEnergyModelVariables{FT,N}
    energybalance::EnergyBalanceVariables{FT,N,Np}
    environmentalconditions::EnvironmentalConditions{FT,N}
    heatflux::HeatFluxVariables{FT,N}
    humidity::HumidityVariables{FT,N}
    radiationflux::RadiationFluxVariables{FT,N}
    temperature::TemperatureVariables{FT,N}
    waterflux::WaterFluxVariables{FT,N,Np}
end

"""
    initialize_model_variable_set(
        ::Type{FT},
        N::Int,
        soil_parameters::SoilParameters{FT},
        vegetation_parameters::VegetationParameters{FT},
        initial_value::FT=400.0,
        hours::Int=1,
    ) where {FT<:AbstractFloat}

Initialize a `ModelVariableSet` with the given parameters.

# Arguments
- `FT::Type`: The floating-point type (e.g., `Float64`).
- `N::Int`: Dimensionality parameter (0 for scalar, 1 for time series).
- `Tatm::FT`: Atmospheric temperature for initialization. Corresponds to `MeteoData.Tatm` in the original MATLAB code
- `AtmSpecific::FT`: Atmospheric specific humidity for initialization. Corresponds to `HumidityAtm.AtmSpecific` in the original MATLAB code.
- `soil_parameters::SoilParameters{FT}`: Soil parameters for initialization.
- `vegetation_parameters::VegetationParameters{FT}`: Vegetation parameters for initialization.
- `initial_value::FT`: Initial value for water flux variables (default is `400.0`).
- `hours::Int`: Number of time steps (default is `1`).
"""
function initialize_model_variable_set(
    ::Type{FT},
    ::TimeSlice,
    soil_parameters::SoilParameters{FT},
    vegetation_parameters::VegetationParameters{FT},
) where {FT<:AbstractFloat}
    N = dimension_value(TimeSlice())
    return initialize(
        FT,
        ModelVariableSet,
        Dict{String,Any}(),
        (FT, N, N+1),
        soil_parameters,
        vegetation_parameters,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{ModelVariableSet},
    data::Dict{String,Any},
    params::Tuple,
    soil_parameters::SoilParameters{FT},
    vegetation_parameters::VegetationParameters{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["buildingenergymodel"] = initialize_building_energy_model_variables(
        FT, dimensionality_type(params[2])
    )
    processed["energybalance"] = initialize_energy_balance_variables(
        FT, dimensionality_type(params[2])
    )
    processed["environmentalconditions"] = initialize_environmental_conditions(
        FT, dimensionality_type(params[2])
    )
    processed["heatflux"] = initialize_heat_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["humidity"] = initialize_humidity_variables(
        FT, dimensionality_type(params[2])
    )
    processed["radiationflux"] = initialize_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["temperature"] = initialize_temperature_variables(
        FT, dimensionality_type(params[2])
    )
    processed["waterflux"] = initialize_water_flux_variables(
        FT, dimensionality_type(params[2]), soil_parameters, vegetation_parameters
    )

    return processed
end

function initialize_model_variable_set(
    ::Type{FT},
    ::TimeSeries,
    Tatm::FT,
    AtmSpecific::FT,
    soil_parameters::SoilParameters{FT},
    vegetation_parameters::VegetationParameters{FT},
    initial_value::FT,
    hours::Int,
) where {FT<:AbstractFloat}
    N = dimension_value(TimeSeries())
    return initialize(
        FT,
        ModelVariableSet,
        Dict{String,Any}(),
        (FT, N, N+1),
        hours,
        Tatm,
        AtmSpecific,
        soil_parameters,
        vegetation_parameters,
        initial_value,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{ModelVariableSet},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    Tatm::FT,
    AtmSpecific::FT,
    soil_parameters::SoilParameters{FT},
    vegetation_parameters::VegetationParameters{FT},
    initial_value::FT,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["buildingenergymodel"] = initialize_building_energy_model_variables(
        FT, dimensionality_type(params[2]), hours, Tatm, AtmSpecific
    )
    processed["energybalance"] = initialize_energy_balance_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["environmentalconditions"] = initialize_environmental_conditions(
        FT, dimensionality_type(params[2]), hours
    )
    processed["heatflux"] = initialize_heat_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["humidity"] = initialize_humidity_variables(
        FT, dimensionality_type(params[2]), hours, AtmSpecific
    )
    processed["radiationflux"] = initialize_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["temperature"] = initialize_temperature_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["waterflux"] = initialize_water_flux_variables(
        FT,
        dimensionality_type(params[2]),
        soil_parameters,
        vegetation_parameters,
        initial_value,
        hours,
    )

    return processed
end
