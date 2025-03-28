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
end

function initialize_vegetated_soilparameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    return initialize(FT, VegetatedSoilParameters, data)
end

function TethysChlorisCore.get_required_fields(::Type{VegetatedSoilParameters})
    return [:Pcla, :Psan, :Porg, :In_max_imp, :Sp_In, :Kimp, :Kfc, :Phy, :SPAR, :Kbot]
end

function TethysChlorisCore.validate_fields(
    ::Type{VegetatedSoilParameters}, data::Dict{String,Any}
)
    # Check that data does not include a key beyond the three components
    for key in keys(data)
        if key ∉ String.(fieldnames(VegetatedSoilParameters))
            throw(ArgumentError("Extraneous key: $key"))
        end
    end
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
    Sp_In_T::FT
end

function initialize_soil_parameters(
    ::Type{FT}, data::Dict{String,Any}
) where {FT<:AbstractFloat}
    processed = copy(data)

    processed["Sp_In_T"] = data["Sp_In_T"]
    processed["roof"] = initialize_vegetated_soilparameters(FT, data["roof"])
    processed["ground"] = initialize_vegetated_soilparameters(FT, data["ground"])

    return initialize(FT, SoilParameters, processed)
end

function TethysChlorisCore.get_required_fields(::Type{SoilParameters})
    return [:roof, :ground, :Sp_In_T]
end

function TethysChlorisCore.validate_fields(::Type{SoilParameters}, data::Dict{String,Any})
    # Check that data does not include a key beyond the four components
    for key in keys(data)
        if key ∉ String.(fieldnames(SoilParameters))
            throw(ArgumentError("Extraneous key: $key"))
        end
    end
end
