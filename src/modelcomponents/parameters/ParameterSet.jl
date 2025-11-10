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
Base.@kwdef struct ParameterSet{FT<:AbstractFloat} <: AbstractParameterSet{FT}
    building_energy::BuildingEnergyModelParameters{FT}
    person::PersonParameters{FT}
    soil::SoilParameters{FT}
    surfacefractions::SurfaceFractions{FT}
    thermal::ThermalProperties{FT}
    optical::OpticalProperties{FT}
    urbangeometry::UrbanGeometryParameters{FT}
    vegetation::VegetationParameters{FT}
    location::LocationProperties{FT}
end

function initialize_parameter_set(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, ParameterSet, data, (FT,))
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{ParameterSet}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    # Initialize each component
    processed["building_energy"] = initialize_building_energy_model_parameters(
        FT, data["building_energy"]
    )
    processed["person"] = initialize_person_parameters(FT, data["person"])

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
    processed["soil"] = initialize_soil_parameters(FT, data["soil"])
    processed["location"] = initialize_locationproperties(FT, data["location"])

    return processed
end

function TethysChlorisCore.validate_fields(::Type{ParameterSet}, data::Dict{String,Any})
    check_extraneous_fields(ParameterSet, data)

    hcanyon, wcanyon, _, htree, radius_tree, distance_tree, _ = preprocess_geometry(
        data["urbangeometry"]["Height_canyon"],
        data["urbangeometry"]["Width_canyon"],
        data["urbangeometry"]["Width_roof"],
        data["urbangeometry"]["Height_tree"],
        data["urbangeometry"]["Radius_tree"],
        data["urbangeometry"]["Distance_tree"],
    )

    # Validate the person position with respect to the urban geometry
    px = data["person"]["PositionPx"] / data["urbangeometry"]["Width_canyon"]
    pz = data["person"]["PositionPz"] / data["urbangeometry"]["Height_canyon"]

    if px >= wcanyon
        throw(ArgumentError("This position of the person P is not permitted"))
    end

    if pz >= hcanyon
        throw(ArgumentError("This position of the person P is not permitted"))
    end

    if data["urbangeometry"]["trees"]
        if sqrt((px - distance_tree)^2 + (pz - radius_tree)^2) <= radius_tree
            throw(ArgumentError("This position of the person P is not permitted"))
        end

        if sqrt(((wcanyon - px) - distance_tree)^2 + (pz - htree)^2) <= radius_tree
            throw(ArgumentError("This position of the person P is not permitted"))
        end
    end
end
