"""
    ModelVariableSet{FT<:AbstractFloat, N, Np} <: Abstract2PModelVariablesSet{FT, N, Np1}

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
Base.@kwdef struct ModelVariableSet{FT<:AbstractFloat,MR,MG} <: AbstractModelVariableSet{FT}
    buildingenergymodel::BuildingEnergyModelVariables{FT}
    energybalance::EnergyBalanceVariables{FT}
    environmentalconditions::EnvironmentalConditions{FT}
    heatflux::HeatFluxVariables{FT}
    humidity::HumidityVariables{FT}
    radiationflux::RadiationFluxVariables{FT}
    temperature::TemperatureVariables{FT}
    waterflux::WaterFluxVariables{FT,MR,MG}
end

"""
    ModelVariableSet(
        ::Type{FT},
        soil_parameters::SoilParameters{FT},
    ) where {FT<:AbstractFloat}

Initialize a `ModelVariableSet` with the given parameters.

# Arguments
- `FT::Type`: The floating-point type (e.g., `Float64`).
- `soil_parameters::SoilParameters{FT}`: Soil parameters for initialization.
"""
function ModelVariableSet(
    ::Type{FT}, soil_parameters::SoilParameters{FT}
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        ModelVariableSet,
        Dict{String,Any}(),
        (FT, soil_parameters.roof.ms, soil_parameters.ground.ms),
        soil_parameters,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{ModelVariableSet},
    data::Dict{String,Any},
    params::Tuple,
    soil_parameters::SoilParameters{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["buildingenergymodel"] = BuildingEnergyModelVariables(FT)
    processed["energybalance"] = EnergyBalanceVariables(FT)
    processed["environmentalconditions"] = EnvironmentalConditions(FT)
    processed["heatflux"] = HeatFluxVariables(FT)
    processed["humidity"] = HumidityVariables(FT)
    processed["radiationflux"] = RadiationFluxVariables(FT)
    processed["temperature"] = TemperatureVariables(FT)
    processed["waterflux"] = WaterFluxVariables(FT, soil_parameters)

    return processed
end

# Automatically define parent accessors for all model variable components
for field in fieldnames(ModelVariableSet)
    component_type = fieldtype(ModelVariableSet, field)
    @eval ModelComponents.parent_accessor(::Type{$component_type}) = x -> x.variables.$field
end

function ModelComponents.accessors(::Type{ModelVariableSet}, ::Type{O}) where {O}
    base = Dict{Symbol,Dict{Symbol,Function}}()

    for field in fieldnames(ModelVariableSet)
        component_type = fieldtype(ModelVariableSet, field)
        component_accessors = ModelComponents.accessors(component_type, O)
        if !isempty(component_accessors)
            merge!(base, component_accessors)
        end
    end
    return base
end

function update!(x::ModelVariableSet, results::NamedTuple, ::EBWBRoofDispatcher)
    update!(x.radiationflux, results, eb_wb_roof_dispatcher)
    update!(x.heatflux, results, eb_wb_roof_dispatcher)
    update!(x.environmentalconditions, results, eb_wb_roof_dispatcher)
    update!(x.waterflux, results, eb_wb_roof_dispatcher)
    update!(x.energybalance, results, eb_wb_roof_dispatcher)

    return nothing
end

function update!(x::ModelVariableSet, results::NamedTuple, ::EBWBCanyonDispatcher)
    update!(x.radiationflux, results, eb_wb_canyon_dispatcher)
    update!(x.heatflux, results, eb_wb_canyon_dispatcher)
    update!(x.environmentalconditions, results, eb_wb_canyon_dispatcher)
    update!(x.waterflux, results, eb_wb_canyon_dispatcher)
    update!(x.energybalance, results, eb_wb_canyon_dispatcher)
    update!(x.humidity, results, eb_wb_canyon_dispatcher)
    update!(x.temperature, results, eb_wb_canyon_dispatcher)

    return nothing
end
