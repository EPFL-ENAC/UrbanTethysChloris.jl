"""
    PersonParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for the person in the urban canyon for MRT calculations.

# Fields
- `PositionPx::FT`: Position within canyon [m]
- `PositionPz::FT`: Height of centre of person [m]
- `PersonWidth::FT`: Horizontal radius of ellipse describing person (=hip width / 2) [-]
- `PersonHeight::FT`: Vertical radius of ellipse describing person (= height / 2) [-]
- `HeightWind::FT`: Height for wind speed to calculate OTC [m]
"""
Base.@kwdef struct PersonParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    PositionPx::FT
    PositionPz::FT
    PersonWidth::FT
    PersonHeight::FT
    HeightWind::FT
end

function initialize_person_parameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, PersonParameters, data)
end

function TethysChlorisCore.get_required_fields(::Type{PersonParameters})
    return [:PositionPx, :PositionPz, :PersonWidth, :PersonHeight, :HeightWind]
end

function TethysChlorisCore.validate_fields(::Type{PersonParameters}, data::Dict{String,Any})
    check_extraneous_fields(
        PersonParameters, data, String.(get_required_fields(PersonParameters))
    )

    # Check if PositionPx and PositionPz are within valid ranges
    if data["PositionPx"] < 0.0 || data["PositionPz"] < 0.0
        throw(ArgumentError("PositionPx and PositionPz must be non-negative"))
    end
end
