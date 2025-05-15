abstract type AbstractAnthropogenicInputs{FT<:AbstractFloat} <: AbstractForcingInputs{FT} end

Base.@kwdef struct AnthropogenicInputs{FT<:AbstractFloat} <: AbstractAnthropogenicInputs{FT}
    Tb::Vector{FT}
    Qf_canyon::Vector{FT}
    Qf_roof::Vector{FT}
    Waterf_canyonVeg::Vector{FT}
    Waterf_canyonBare::Vector{FT}
    Waterf_roof::Vector{FT}
end

function initialize_anthropogenic_inputs(
    ::Type{FT}, data::NCDataset, Tatm::Vector{FT}
) where {FT<:AbstractFloat}
    return initialize(FT, AnthropogenicInputs, data, Tatm)
end

function TethysChlorisCore.get_required_fields(::Type{AnthropogenicInputs})
    return [:Tbmin, :Tbmax]
end

function TethysChlorisCore.get_optional_fields(::Type{AnthropogenicInputs})
    return [:Qf_canyon, :Qf_roof, :Waterf_canyonVeg, :Waterf_canyonBare, :Waterf_roof]
end

function TethysChlorisCore.get_calculated_fields(::Type{AnthropogenicInputs})
    return [:Tb]
end

function get_array_fields(::Type{AnthropogenicInputs})
    return [:Tb, :Qf_canyon, :Qf_roof, :Waterf_canyonVeg, :Waterf_canyonBare, :Waterf_roof]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{AnthropogenicInputs}, data::NCDataset, Tatm::Vector{FT}
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    for field in TethysChlorisCore.get_optional_fields(AnthropogenicInputs)
        if haskey(data, field)
            if dimnames(data[field]) == ()
                processed[String(field)] = fill(FT(data[field][]), data.dim["hours"])
            else
                processed[String(field)] = Array(data[field])
            end
        else
            processed[String(field)] = zeros(FT, data.dim["hours"])
        end
    end

    # Calculate Tb
    processed["Tb"] = Tb(Tatm, data["Tbmin"][], data["Tbmax"][])

    return processed
end

function Tb(Tatm::Vector{FT}, Tbmin::FT, Tbmax::FT) where {FT<:AbstractFloat}
    return clamp.(Tatm, Tbmin + 273.15, Tbmax + 273.15)
end

function TethysChlorisCore.validate_fields(::Type{AnthropogenicInputs}, data::NCDataset) end
