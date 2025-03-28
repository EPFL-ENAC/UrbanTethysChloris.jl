abstract type AbstractParameterSet{FT<:AbstractFloat} <: AbstractModelComponentSet{FT} end
const APS = AbstractParameterSet

"""
    ParameterSet{FT<:AbstractFloat} <: AbstractParameterSet{FT}

Parameters for the Urban Tethys-Chloris model.

# Fields
- `building_energy::BuildingEnergyModelParameters{FT}`: Parameters for the building energy model.
- `person::PersonParameters{FT}`: Parameters for the person.
- `soil::SoilParameters{FT}`: Parameters for the soil.
- `surfacefractions::SurfaceFractions{FT}`: Parameters for the surface fractions.
- `thermal::ThermalProperties{FT}`: Parameters for the thermal properties.
- `optical::OpticalProperties{FT}`: Parameters for the optical properties.
- `urbangeometry::UrbanGeometryParameters{FT}`: Parameters for the urban geometry.
- `vegetation::VegetationParameters{FT}`: Parameters for the vegetation.
"""
Base.@kwdef struct ParameterSet{FT<:AbstractFloat} <: APS{FT}
    building_energy::BuildingEnergyModelParameters{FT}
    person::PersonParameters{FT}
    soil::SoilParameters{FT}
    surfacefractions::SurfaceFractions{FT}
    thermal::ThermalProperties{FT}
    optical::OpticalProperties{FT}
    urbangeometry::UrbanGeometryParameters{FT}
    vegetation::VegetationParameters{FT}
end

function initialize_parameter_set(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, ParameterSet, data)
end

function TethysChlorisCore.get_required_fields(::Type{ParameterSet})
    return [
        :building_energy,
        :person,
        :soil,
        :surfacefractions,
        :thermal,
        :optical,
        :urbangeometry,
        :vegetation,
    ]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{ParameterSet}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Initialize each component
    processed["building_energy"] = initialize_building_energy_model_parameters(
        FT, data["building_energy"]
    )
    processed["person"] = initialize_person_parameters(FT, data["person"])
    processed["soil"] = initialize_soil_parameters(FT, data["soil"])

    processed["surfacefractions"] = initialize_surfacefractions(
        FT, data["surfacefractions"]
    )
    processed["thermal"] = initialize_thermalproperties(FT, data["thermal"])
    processed["optical"] = initialize_optical_properties(
        FT, data["optical"], processed["surfacefractions"]
    )
    processed["urbangeometry"] = initialize_urbangeometry_parameters(
        FT, data["urbangeometry"]
    )
    processed["vegetation"] = initialize_vegetationparameters(FT, data["vegetation"])

    return processed
end
