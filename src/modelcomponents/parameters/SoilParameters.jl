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
    Zs::Vector{FT}
    ms::Int
    FixSM::Bool
    FixSM_LayerStart::Int
    FixSM_LayerEnd::Int
    Ohy::FT
    Osat::FT
    O33::FT
    dz::Vector{FT}
end

function initialize_vegetated_soilparameters(
    ::Type{FT},
    data::Dict{String,Any},
    tree::HeightDependentVegetationParameters{FT},
    ground::HeightDependentVegetationParameters{FT},
) where {FT<:AbstractFloat}
    return initialize(FT, VegetatedSoilParameters, data, (FT,), tree, ground)
end

function TethysChlorisCore.get_calculated_fields(::Type{VegetatedSoilParameters})
    return [:ms, :Osat, :Ohy, :O33, :dz]
end

function TethysChlorisCore.get_optional_fields(::Type{VegetatedSoilParameters})
    return [:In_max_ground, :In_max_underveg, :In_max_bare, :dz1, :dz2]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{VegetatedSoilParameters},
    data::Dict{String,Any},
    params::Tuple,
    tree::HeightDependentVegetationParameters{FT},
    ground::HeightDependentVegetationParameters{FT},
) where {FT<:AbstractFloat}
    processed = copy(data)

    processed["ms"] = length(processed["Zs"]) - 1

    _, _, _, Osat, Ohy, _, _, _, _, _, O33 = Soil.soil_parameters_total(
        processed["Pcla"],
        processed["Psan"],
        processed["Porg"],
        processed["Kfc"],
        processed["Phy"],
        processed["SPAR"],
        processed["Kbot"],
        tree.CASE_ROOT,
        ground.CASE_ROOT,
        tree.ZR95,
        ground.ZR95,
        tree.ZR50,
        ground.ZR50,
        tree.ZRmax,
        ground.ZRmax,
        processed["Zs"],
    )

    unique_Osat = unique(Osat)
    @assert length(unique_Osat) == 1 "Osat should be unique after preprocessing"
    unique_Ohy = unique(Ohy)
    @assert length(unique_Ohy) == 1 "Ohy should be unique after preprocessing"
    unique_O33 = unique(O33)
    @assert length(unique_O33) == 1 "O33 should be unique after preprocessing"

    processed["Osat"] = unique_Osat[1]
    processed["Ohy"] = unique_Ohy[1]
    processed["O33"] = unique_O33[1]
    processed["dz"] = diff(processed["Zs"])

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
    ::Type{FT}, data::Dict{String,Any}, vegetation::VegetationParameters{FT}
) where {FT<:AbstractFloat}
    processed = copy(data)

    processed["roof"] = initialize_vegetated_soilparameters(
        FT, data["roof"], vegetation.roof, vegetation.roof
    )
    processed["ground"] = initialize_vegetated_soilparameters(
        FT, data["ground"], vegetation.tree, vegetation.ground
    )
    processed["wall"] = initialize_wall_soilparameters(FT, data["wall"])

    return initialize(FT, SoilParameters, processed)
end
