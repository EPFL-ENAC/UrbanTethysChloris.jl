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
