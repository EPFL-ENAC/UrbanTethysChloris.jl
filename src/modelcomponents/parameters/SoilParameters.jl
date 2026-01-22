"""
    VegetatedSoilParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Soil parameters specific to vegetated surfaces (roof and ground)

# Fields
- `Pcla::FT`: Fraction of clay in the soil (-)
- `Psan::FT`: Fraction of sand in the soil (-)
- `Porg::FT`: Fraction of organic material in the soil (-)
- `In_max_imp::FT`: Maxiumum interception capacity of impervious area (mm)
- `In_max_ground::FT`: Maxiumum interception capacity of ground under roof vegetation (mm)
- `In_max_underveg::FT`: Maxiumum interception capacity of vegetated ground area (mm)
- `In_max_bare::FT`: Maxiumum interception capacity of bare ground area (mm)
- `Sp_In::FT`: specific water retained by a vegetated surface (mm m^2 VEG area m^-2 plant area)
- `Kimp::FT`: Hydraulic conductivity of impervious area (mm/h)
- `Kfc::FT`: Conductivity at field capacity (mm/h)
- `Phy::FT`: Suction at the residual/hygroscopic water content (kPa)
- `SPAR::Int`: SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
- `Kbot::FT`: Conductivity at the bedrock layer (mm/h)
"""
Base.@kwdef struct VegetatedSoilParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    Pcla::FT
    Psan::FT
    Porg::FT
    In_max_imp::FT
    In_max_ground::FT = FT(NaN)
    In_max_underveg::FT = FT(NaN)
    In_max_bare::FT = FT(NaN)
    Sp_In::FT
    Kimp::FT
    Kfc::FT
    Phy::FT
    SPAR::Int
    Kbot::FT
    dz1::FT = FT(NaN)
    dz2::FT = FT(NaN)
    Zs::Vector{FT} # TODO: replace with static array
    ms::Int
    FixSM::Bool
    FixSM_LayerStart::Int
    FixSM_LayerEnd::Int
end

function initialize_vegetated_soilparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, VegetatedSoilParameters, data, (FT,))
end

function VegetatedSoilParameters(
    ::Type{FT}, data::AbstractDict, suffix::String=""
) where {FT<:AbstractFloat}
    # remove suffix from keys
    stripped_data = Dict{String,Any}()
    for (k, v) in data
        if endswith(k, suffix)
            new_key = replace(k, suffix => "")
            stripped_data[new_key] = v
        end
    end
    return VegetatedSoilParameters{FT}(
        data["Pcla"],
        data["Psan"],
        data["Porg"],
        data["In_max_imp"],
        get(data, "In_max_ground", FT(NaN)),
        get(data, "In_max_underveg", FT(NaN)),
        get(data, "In_max_bare", FT(NaN)),
        data["Sp_In"],
        data["Kimp"],
        data["Kfc"],
        data["Phy"],
        data["SPAR"],
        data["Kbot"],
        get(data, "dz1", FT(NaN)),
        get(data, "dz2", FT(NaN)),
        data["Zs"],
        data["ms"],
        stripped_data["FixSM"],
        stripped_data["FixSM_LayerStart"],
        stripped_data["FixSM_LayerEnd"],
    )
end

function TethysChlorisCore.get_calculated_fields(::Type{VegetatedSoilParameters})
    return [:ms]
end

function TethysChlorisCore.get_optional_fields(::Type{VegetatedSoilParameters})
    return [:In_max_ground, :In_max_underveg, :In_max_bare, :dz1, :dz2]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{VegetatedSoilParameters}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = copy(data)

    processed["ms"] = length(processed["Zs"]) - 1

    return processed
end

Base.@kwdef struct WallSoilParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    dz1::FT
    dz2::FT
end

function initialize_wall_soilparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, WallSoilParameters, data)
end

"""
    SoilParameters{FT<:AbstractFloat} <: AbstractParameters{FT}

Container for soil parameters for different urban surface components.

# Fields
- `roof::VegetatedSoilParameters{FT}`: Roof soil parameters
- `ground::VegetatedSoilParameters{FT}`: Ground soil parameters
- `Sp_In_T::FT`: Specific water retained by a tree (mm m^2 VEG area m^-2 plant area)
"""
Base.@kwdef struct SoilParameters{FT<:AbstractFloat} <: AbstractParameters{FT}
    roof::VegetatedSoilParameters{FT}
    ground::VegetatedSoilParameters{FT}
    wall::WallSoilParameters{FT}
    Sp_In_T::FT
end

function initialize_soil_parameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = copy(data)

    processed["roof"] = initialize_vegetated_soilparameters(FT, data["roof"])
    processed["ground"] = initialize_vegetated_soilparameters(FT, data["ground"])
    processed["wall"] = initialize_wall_soilparameters(FT, data["wall"])

    return initialize(FT, SoilParameters, processed)
end
