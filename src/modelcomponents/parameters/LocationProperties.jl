Base.@kwdef struct LocationProperties{FT<:AbstractFloat} <: AbstractParameters{FT}
    phi::FT
    lambda::FT
    theta_canyon::FT
    DeltaGMT::FT
end

function initialize_locationproperties(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, LocationProperties, data, (FT,))
end

function TethysChlorisCore.get_required_fields(::Type{LocationProperties})
    return [:phi, :lambda, :theta_canyon, :DeltaGMT]
end
