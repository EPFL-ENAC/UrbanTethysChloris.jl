"""
    UrbanGeometryParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Parameters for the UrbanGeometry model component.

# Fields
- `fveg_R::FT`: Vegetated roof fraction [-]
- `fimp_R::FT`: Impervious roof fraction [-].
- `Per_runoff_R::FT`: Fraction of excess water that leaves the system as runoff [-].
- `fveg_G::FT`: Vegetated ground fraction [-].
- `fbare_G::FT`: Bare ground fraction [-].
- `fimp_G::FT`: Impervious ground fraction [-].
- `Per_runoff_G::FT`: Fraction of excess water that leaves the system as runoff [-].
"""
Base.@kwdef struct SurfaceFractions{FT<:AbstractFloat} <: AbstractParameters{FT}
    fveg_R::FT
    fimp_R::FT
    Per_runoff_R::FT
    fveg_G::FT
    fbare_G::FT
    fimp_G::FT
    Per_runoff_G::FT
end

function initialize_surfacefractions(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, SurfaceFractions, data)
end

function TethysChlorisCore.get_required_fields(::Type{SurfaceFractions})
    return [:fveg_R, :fimp_R, :Per_runoff_R, :fveg_G, :fbare_G, :fimp_G, :Per_runoff_G]
end

function TethysChlorisCore.validate_fields(::Type{SurfaceFractions}, data::Dict{String,Any})
    if data["fveg_R"] + data["fimp_R"] != 1.0
        throw(ArgumentError("Roof fractions must sum to 1.0"))
    end

    if data["fveg_G"] + data["fbare_G"] + data["fimp_G"] != 1.0
        throw(ArgumentError("Ground fractions must sum to 1.0"))
    end

    if data["Per_runoff_R"] < 0.0 || data["Per_runoff_R"] > 1.0
        throw(ArgumentError("Roof runoff fraction must be between 0.0 and 1.0"))
    end

    if data["Per_runoff_G"] < 0.0 || data["Per_runoff_G"] > 1.0
        throw(ArgumentError("Ground runoff fraction must be between 0.0 and 1.0"))
    end

    # return nothing
end
