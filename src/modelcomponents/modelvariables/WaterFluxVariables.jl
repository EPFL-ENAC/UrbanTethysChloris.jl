"""
    Eflux{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Evapotranspiration fluxes.

# Fields
- `EfluxRoofImp`: Evaporation flux of intercepted water from impervious roof area to atmosphere [kg/m²*s horizontal impervious roof area]
- `EfluxRoofVegInt`: Evaporation flux of intercepted water on roof vegetation to atmosphere [kg/m²*s horizontal vegetated roof area]
- `EfluxRoofVegPond`: Evaporation flux of intercepted water on ground under roof vegetation to atmosphere [kg/m²*s horizontal vegetated roof area]
- `EfluxRoofVegSoil`: Evaporation flux of water from roof soil under vegetation to atmosphere [kg/m²*s horizontal vegetated roof area]
- `TEfluxRoofVeg`: Evaporation flux of transpiration from roof plants to atmosphere [kg/m²*s horizontal vegetated roof area]
- `EfluxRoofVeg`: Total evaporation flux of vegetated roof to atmosphere [kg/m²*s horizontal vegetated roof area]
- `EfluxRoof`: Total evaporation flux of roof to atmosphere [kg/m²*s horizontal vegetated roof area]
- `EfluxGroundImp`: Evaporation flux of intercepted water on impervious ground area to canyon [kg/m²*s horizontal impervious ground area]
- `EfluxGroundBarePond`: Evaporation flux of water on bare ground to canyon [kg/m²*s horizontal bare ground area]
- `EfluxGroundBareSoil`: Evaporation flux of water from bare ground to canyon [kg/m²*s horizontal bare ground area]
- `EfluxGroundBare`: Total evaporation flux of bare ground area to canyon [kg/m²*s horizontal bare ground area]
- `EfluxGroundVegInt`: Evaporation flux of intercepted water on ground vegetation to canyon [kg/m²*s horizontal vegetated ground area]
- `EfluxGroundVegPond`: Evaporation flux of intercepted water on ground under vegetation to canyon [kg/m²*s horizontal vegetated ground area]
- `EfluxGroundVegSoil`: Evaporation flux of water from ground soil under vegetation to canyon [kg/m²*s horizontal vegetated ground area]
- `TEfluxGroundVeg`: Evaporation flux of transpiration from ground plants to canyon [kg/m²*s horizontal vegetated ground area]
- `EfluxGroundVeg`: Total evaporation flux of vegetated ground to canyon [kg/m²*s horizontal vegetated ground area]
- `EfluxGround`: Total evaporation flux of ground to canyon [kg/m²*s horizontal vegetated ground area]
- `EfluxTreeInt`: Evaporation flux of intercepted water on tree canopy to canyon [kg/m²*s horizontally projected tree area]
- `TEfluxTree`: Evaporation flux of transpiration from tree canopy to canyon [kg/m²*s horizontally projected tree area]
- `EfluxTree`: Total evaporation flux of tree canopy to canyon [kg/m²*s horizontally projected tree area]
- `EfluxWallSun`: Evaporation flux of sunlit wall to canyon [kg/m²*s vertical wall area]
- `EfluxWallShade`: Evaporation flux of shaded wall to canyon [kg/m²*s vertical wall area]
- `EfluxCanyon`: Evaporation flux of canyon to atmosphere [kg/m²*s horizontal canyon area]
- `EfluxUrban`: Total evaporation flux of urban area to atmosphere [kg/m²*s horizontal urban area]
"""
Base.@kwdef mutable struct Eflux{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    EfluxRoofImp::Array{FT,N}
    EfluxRoofVegInt::Array{FT,N}
    EfluxRoofVegPond::Array{FT,N}
    EfluxRoofVegSoil::Array{FT,N}
    TEfluxRoofVeg::Array{FT,N}
    EfluxRoofVeg::Array{FT,N}
    EfluxRoof::Array{FT,N}
    EfluxGroundImp::Array{FT,N}
    EfluxGroundBarePond::Array{FT,N}
    EfluxGroundBareSoil::Array{FT,N}
    EfluxGroundBare::Array{FT,N}
    EfluxGroundVegInt::Array{FT,N}
    EfluxGroundVegPond::Array{FT,N}
    EfluxGroundVegSoil::Array{FT,N}
    TEfluxGroundVeg::Array{FT,N}
    EfluxGroundVeg::Array{FT,N}
    EfluxGround::Array{FT,N}
    EfluxTreeInt::Array{FT,N}
    TEfluxTree::Array{FT,N}
    EfluxTree::Array{FT,N}
    EfluxWallSun::Array{FT,N}
    EfluxWallShade::Array{FT,N}
    EfluxCanyon::Array{FT,N}
    EfluxUrban::Array{FT,N}
end

function initialize_eflux(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Eflux, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_eflux(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, Eflux, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    Runoff{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Runoff variables for different urban surfaces.

# Fields
- `QRoofImp`: Runoff of impervious area of roof [mm/time step per horizontal impervious roof area]
- `QRoofVegDrip`: Runoff, Dripping, etc from vegetation to roof ground [mm/time step per horizontal vegetated roof area]
- `QRoofVegPond`: Runoff from roof soil under vegetation due to limitation in infiltration capacity [mm/time step per horizontal vegetated roof area]
- `QRoofVegSoil`: Runoff due to roof soil saturation [mm/time step per horizontal vegetated roof area]
- `QGroundImp`: Runoff of impervious area of ground [mm/time step per horizontal impervious ground area]
- `QGroundBarePond`: Runoff of bare area of ground due to limitation in infiltration capacity [mm/time step per horizontal bare ground area]
- `QGroundBareSoil`: Runoff of bare area of ground due to soil saturation [mm/time step per horizontal bare ground area]
- `QTree`: Runoff, Dripping, etc from tree to ground [mm/time step per horizontally projected tree area]
- `QGroundVegDrip`: Runoff, Dripping, etc from vegetation to ground [mm/time step per horizontal vegetated ground area]
- `QGroundVegPond`: Runoff from soil under vegetation due to limitation in infiltration capacity [mm/time step per horizontal vegetated ground area]
- `QGroundVegSoil`: Runoff due to soil saturation under vegetation on ground [mm/time step per horizontal vegetated ground area]
"""
Base.@kwdef mutable struct Runoff{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    QRoofImp::Array{FT,N}
    QRoofVegDrip::Array{FT,N}
    QRoofVegPond::Array{FT,N}
    QRoofVegSoil::Array{FT,N}
    QGroundImp::Array{FT,N}
    QGroundBarePond::Array{FT,N}
    QGroundBareSoil::Array{FT,N}
    QTree::Array{FT,N}
    QGroundVegDrip::Array{FT,N}
    QGroundVegPond::Array{FT,N}
    QGroundVegSoil::Array{FT,N}
end

function initialize_runoff(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Runoff, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_runoff(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, Runoff, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    Runon{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Runon variables for urban area.

# Fields
- `RunonRoofTot`: Total roof runon to the next time step [mm/time step per horizontal roof area]
- `RunoffRoofTot`: Total roof runoff that is removed from the system [mm/time step per horizontal roof area]
- `RunonGroundTot`: Total runon in canyon to the next time step [mm/time step per horizontal ground area]
- `RunoffGroundTot`: Total runoff in canyon that is removed from the system [mm/time step per horizontal ground area]
- `RunonUrban`: Total urban runon to the next time step [mm/time step per horizontal urban area]
- `RunoffUrban`: Total urban runoff that is removed from the system [mm/time step per horizontal urban area]
"""
Base.@kwdef mutable struct Runon{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    RunonRoofTot::Array{FT,N}
    RunoffRoofTot::Array{FT,N}
    RunonGroundTot::Array{FT,N}
    RunoffGroundTot::Array{FT,N}
    RunonUrban::Array{FT,N}
    RunoffUrban::Array{FT,N}
end

function initialize_runon(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Runon, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_runon(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, Runon, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    Leakage{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Leakage variables for different urban surfaces.

# Fields
- `LkRoofImp`: Leakage from impervious roof [mm/h per horizontal impervious roof area]
- `LkRoofVeg`: Leakage from last soil layer of vegetated roof [mm/h per horizontal vegetated roof area]
- `LkRoof`: Total leakage of roof [mm/h per horizontal roof area]
- `LkGroundImp`: Leakage from impervious ground [mm/h per horizontal impervious ground area]
- `LkGroundBare`: Leakage from last soil layer of bare ground [mm/h per horizontal bare ground area]
- `LkGroundVeg`: Leakage from last soil layer of vegetated ground [mm/h per horizontal vegetated ground area]
- `LkGround`: Total leakage of ground [mm/h per horizontal ground area]
- `LkUrban`: Total leakage of ground and roof soil [mm/h per horizontal urban area]
"""
Base.@kwdef mutable struct Leakage{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    LkRoofImp::Array{FT,N}
    LkRoofVeg::Array{FT,N}
    LkRoof::Array{FT,N}
    LkGroundImp::Array{FT,N}
    LkGroundBare::Array{FT,N}
    LkGroundVeg::Array{FT,N}
    LkGround::Array{FT,N}
    LkUrban::Array{FT,N}
end

function initialize_leakage(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Leakage, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_leakage(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, Leakage, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    Interception{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Interception variables for different urban surfaces.

# Fields
- `IntRoofImp`: Interception on impervious roof area [mm per horizontal impervious roof area]
- `IntRoofVegPlant`: Interception on plant surfaces [mm per horizontal vegetated roof area]
- `IntRoofVegGround`: Interception on ground [mm per horizontal roof area]
- `IntRooftot`: Total interception on roof [mm per horizontal roof area]
- `IntGroundImp`: Interception on impervious ground area [mm per horizontal impervious ground area]
- `IntGroundBare`: Interception on bare ground area [mm per horizontal bare ground area]
- `IntGroundVegPlant`: Interception on plant surfaces [mm per horizontal vegetated ground area]
- `IntGroundVegGround`: Interception on ground [mm per horizontal vegetated ground area]
- `IntTree`: Interception on tree [mm per horizontally projected tree area]
"""
Base.@kwdef mutable struct Interception{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    IntRoofImp::Array{FT,N}
    IntRoofVegPlant::Array{FT,N}
    IntRoofVegGround::Array{FT,N}
    IntRooftot::Array{FT,N}
    IntGroundImp::Array{FT,N}
    IntGroundBare::Array{FT,N}
    IntGroundVegPlant::Array{FT,N}
    IntGroundVegGround::Array{FT,N}
    IntTree::Array{FT,N}
end

function initialize_interception(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, Interception, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_interception(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT, Interception, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    dInt_dt{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Change in interception variables for different urban surfaces.

# Fields
- `dInt_dtRoofImp`: Change in interception on impervious roof area [mm/h per horizontal impervious roof area]
- `dInt_dtRoofVegPlant`: Change in interception on plant surfaces [mm/h per horizontal vegetated roof area]
- `dInt_dtRoofVegGround`: Change in interception on ground [mm/h per horizontal vegetated roof area]
- `dInt_dtRooftot`: Total change in interception on roof [mm/h per horizontal roof area]
- `dInt_dtGroundImp`: Change in interception on impervious ground area [mm/h per horizontal impervious ground area]
- `dInt_dtGroundBare`: Change in interception on bare ground area [mm/h per horizontal bare ground area]
- `dInt_dtGroundVegPlant`: Change in interception on plant surfaces [mm/h per horizontal vegetated ground area]
- `dInt_dtGroundVegGround`: Change in interception on ground [mm/h per horizontal vegetated ground area]
- `dInt_dtTree`: Change in interception on tree [mm/h per horizontally projected tree area]
"""
# Same names as Interception but with dInt_dt prefix
Base.@kwdef mutable struct dInt_dt{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    dInt_dtRoofImp::Array{FT,N}
    dInt_dtRoofVegPlant::Array{FT,N}
    dInt_dtRoofVegGround::Array{FT,N}
    dInt_dtRooftot::Array{FT,N}
    dInt_dtGroundImp::Array{FT,N}
    dInt_dtGroundBare::Array{FT,N}
    dInt_dtGroundVegPlant::Array{FT,N}
    dInt_dtGroundVegGround::Array{FT,N}
    dInt_dtTree::Array{FT,N}
end

function initialize_dint_dt(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, dInt_dt, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_dint_dt(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, dInt_dt, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    Infiltration{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Infiltration variables for different urban surfaces.

# Fields
- `fRoofVeg`: Infiltration in first soil layer of vegetated roof [mm/h per horizontal roof area]
- `fGroundBare`: Infiltration in first soil layer of bare ground [mm/h per horizontal bare ground area]
- `fGroundVeg`: Infiltration in first soil layer of vegetated ground [mm/h per horizontal vegetated ground area]
- `fGroundImp`: Infiltration in impervious ground (usually zero) [mm/h per horizontal impervious ground area]
"""
Base.@kwdef mutable struct Infiltration{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    fRoofVeg::Array{FT,N}
    fGroundBare::Array{FT,N}
    fGroundVeg::Array{FT,N}
    fGroundImp::Array{FT,N}
end

function initialize_infiltration(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, Infiltration, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_infiltration(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT, Infiltration, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

abstract type AbstractLayeredSoilVariables{FT<:AbstractFloat,N} <:
              Abstract1PModelVariables{FT,N} end

function allocate_fields!(
    x::Dict{String,Any},
    ::Type{FT},
    ::Type{T},
    fields::Vector{String},
    vector_length::Signed,
) where {FT<:AbstractFloat,T<:AbstractLayeredSoilVariables}
    for field in fields
        x[field] = fill(zero(MVector{vector_length,FT}))
    end
end

function allocate_fields!(
    x::Dict{String,Any},
    ::Type{FT},
    ::Type{T},
    fields::Vector{String},
    vector_length::Signed,
    hours::Signed,
) where {FT<:AbstractFloat,T<:AbstractLayeredSoilVariables}
    for field in fields
        x[field] = [zeros(MVector{vector_length,FT}) for _ in 1:hours]
    end
end

# Default initialization does nothing
function initialize_fields!(
    x::Dict{String,Any}, ::Type{T}, fields::Vector{String}, soil_values::NamedTuple, args...
) where {T<:AbstractLayeredSoilVariables} end

function add_fields!(
    x::Dict{String,Any},
    ::Type{FT},
    ::Type{T},
    fields::Vector{String},
    soil_values::NamedTuple,
    args...,
) where {FT<:AbstractFloat,T<:AbstractLayeredSoilVariables}
    allocate_fields!(x, FT, T, fields, soil_values.ms, args...)

    initialize_fields!(x, T, fields, soil_values, args...)
end

function roof_fields(::Type{T}) where {T<:AbstractLayeredSoilVariables}
    return String[]
end

function ground_fields(::Type{T}) where {T<:AbstractLayeredSoilVariables}
    return String[]
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    soil_values::NamedTuple,
    args...,
) where {FT<:AbstractFloat,T<:AbstractLayeredSoilVariables}
    processed = Dict{String,Any}()

    add_fields!(processed, FT, T, ground_fields(T), soil_values.ground, args...)

    add_fields!(processed, FT, T, roof_fields(T), soil_values.roof, args...)

    return processed
end

"""
    Vwater{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Water volume in soil for different urban surfaces.

# Fields
- `VRoofSoilVeg`: Water volume in the different soil layers of roof [mm per horizontal roof area]
- `VGroundSoilImp`: Water volume in the different soil layers of ground under impervious [mm per horizontal impervious ground area]
- `VGroundSoilBare`: Water volume in the different soil layers of ground under bare [mm per horizontal bare ground area]
- `VGroundSoilVeg`: Water volume in the different soil layers of ground under vegetated [mm per horizontal vegetated ground area]
- `VGroundSoilTot`: Water volume in the different soil layers of ground total [mm per horizontal ground area]
"""
Base.@kwdef mutable struct Vwater{FT<:AbstractFloat,N,MR,MG} <:
                           AbstractLayeredSoilVariables{FT,N}
    VRoofSoilVeg::Array{MVector{MR,FT},N}
    VGroundSoilImp::Array{MVector{MG,FT},N}
    VGroundSoilBare::Array{MVector{MG,FT},N}
    VGroundSoilVeg::Array{MVector{MG,FT},N}
    VGroundSoilTot::Array{MVector{MG,FT},N}
end

function initialize_vwater(
    ::Type{FT}, dim::T, soil_values::NamedTuple, args...
) where {FT<:AbstractFloat,T<:ModelDimension}
    return initialize(
        FT,
        Vwater,
        Dict{String,Any}(),
        (FT, dimension_value(dim), soil_values.roof.ms, soil_values.ground.ms),
        soil_values,
        args...,
    )
end

function ground_fields(::Type{Vwater})
    return ["VGroundSoilImp", "VGroundSoilBare", "VGroundSoilVeg", "VGroundSoilTot"]
end

function roof_fields(::Type{Vwater})
    return ["VRoofSoilVeg"]
end

function initialize_fields!(
    x::Dict{String,Any},
    ::Type{Vwater},
    fields::Vector{String},
    soil_values::NamedTuple,
    hours::Signed,
)
    inital_values = soil_values.O33 * soil_values.dz

    for field in fields
        x[field][1][:] = inital_values[:]
    end
end

"""
    dVwater_dt{FT<:AbstractFloat, N, MR, MG} <: AbstractLayeredSoilVariables{FT,N}

Change in water volume in soil for different urban surfaces.

# Fields
- `dVRoofSoilVeg_dt`: Change in water volume in the different soil layers of roof [mm per horizontal roof area]
- `dVGroundSoilImp_dt`: Change in water volume in the different soil layers of ground under impervious [mm per horizontal impervious ground area]
- `dVGroundSoilBare_dt`: Change in water volume in the different soil layers of ground under bare [mm per horizontal bare ground area]
- `dVGroundSoilVeg_dt`: Change in water volume in the different soil layers of ground under vegetated [mm per horizontal ground area]
- `dVGroundSoilTot_dt`: Change in water volume in the different soil layers of ground total [mm per horizontal ground area]
"""
# Same names as Vwater
Base.@kwdef mutable struct dVwater_dt{FT<:AbstractFloat,N,MR,MG} <:
                           AbstractLayeredSoilVariables{FT,N}
    dVRoofSoilVeg_dt::Array{MVector{MR,FT},N}
    dVGroundSoilImp_dt::Array{MVector{MG,FT},N}
    dVGroundSoilBare_dt::Array{MVector{MG,FT},N}
    dVGroundSoilVeg_dt::Array{MVector{MG,FT},N}
    dVGroundSoilTot_dt::Array{MVector{MG,FT},N}
end

function ground_fields(::Type{dVwater_dt})
    return [
        "dVGroundSoilImp_dt",
        "dVGroundSoilBare_dt",
        "dVGroundSoilVeg_dt",
        "dVGroundSoilTot_dt",
    ]
end

function roof_fields(::Type{dVwater_dt})
    return ["dVRoofSoilVeg_dt"]
end

function initialize_dvwater_dt(
    ::Type{FT}, dim::T, soil_values::NamedTuple, args...
) where {FT<:AbstractFloat,T<:ModelDimension}
    return initialize(
        FT,
        dVwater_dt,
        Dict{String,Any}(),
        (FT, dimension_value(dim), soil_values.roof.ms, soil_values.ground.ms),
        soil_values,
        args...,
    )
end

"""
    Owater{FT<:AbstractFloat, N, MR, MG} <: AbstractLayeredSoilVariables{FT,N}

Soil moisture in different soil layers for urban surfaces.

# Fields
- `OwRoofSoilVeg`: Soil moisture in the different soil layers of roof [-]
- `OwGroundSoilImp`: Soil moisture in the different soil layers of ground under impervious [-]
- `OwGroundSoilBare`: Soil moisture in the different soil layers of ground under bare [-]
- `OwGroundSoilVeg`: Soil moisture in the different soil layers of ground under vegetated [-]
- `OwGroundSoilTot`: Soil moisture in the different soil layers of ground total [-]
"""
# Same names as Vwater
Base.@kwdef mutable struct Owater{FT<:AbstractFloat,N,MR,MG} <:
                           AbstractLayeredSoilVariables{FT,N}
    OwRoofSoilVeg::Array{MVector{MR,FT},N}
    OwGroundSoilImp::Array{MVector{MG,FT},N}
    OwGroundSoilBare::Array{MVector{MG,FT},N}
    OwGroundSoilVeg::Array{MVector{MG,FT},N}
    OwGroundSoilTot::Array{MVector{MG,FT},N}
end

function roof_fields(::Type{Owater})
    return ["OwRoofSoilVeg"]
end

function ground_fields(::Type{Owater})
    return ["OwGroundSoilImp", "OwGroundSoilBare", "OwGroundSoilVeg", "OwGroundSoilTot"]
end

function initialize_fields!(
    x::Dict{String,Any},
    ::Type{Owater},
    fields::Vector{String},
    soil_values::NamedTuple,
    hours::Signed,
)
    for field in fields
        x[field][1][:] .= soil_values.O33
    end
end

function initialize_owater(
    ::Type{FT}, dim::T, soil_values::NamedTuple, args...
) where {FT<:AbstractFloat,T<:ModelDimension}
    return initialize(
        FT,
        Owater,
        Dict{String,Any}(),
        (FT, dimension_value(dim), soil_values.roof.ms, soil_values.ground.ms),
        soil_values,
        args...,
    )
end

"""
    OSwater{FT<:AbstractFloat, N} <: AbstractLayeredSoilVariables{FT,N}

Additional soil moisture variables for urban surfaces.

# Fields
- `OSwRoofSoilVeg`: Additional soil moisture values for roof soil layers [-]
- `OSwGroundSoilImp`: Additional soil moisture values for ground soil layers under impervious [-]
- `OSwGroundSoilBare`: Additional soil moisture values for ground soil layers under bare [-]
- `OSwGroundSoilVeg`: Additional soil moisture values for ground soil layers under vegetated [-]
- `OSwGroundSoilTot`: Additional soil moisture values for ground soil layers total [-]
"""
# Same names as Vwater
Base.@kwdef mutable struct OSwater{FT<:AbstractFloat,N,MR,MG} <:
                           AbstractLayeredSoilVariables{FT,N}
    OSwRoofSoilVeg::Array{MVector{MR,FT},N}
    OSwGroundSoilImp::Array{MVector{MG,FT},N}
    OSwGroundSoilBare::Array{MVector{MG,FT},N}
    OSwGroundSoilVeg::Array{MVector{MG,FT},N}
    OSwGroundSoilTot::Array{MVector{MG,FT},N}
end

function roof_fields(::Type{OSwater})
    return ["OSwRoofSoilVeg"]
end

function ground_fields(::Type{OSwater})
    return ["OSwGroundSoilImp", "OSwGroundSoilBare", "OSwGroundSoilVeg", "OSwGroundSoilTot"]
end

function initialize_oswater(
    ::Type{FT}, dim::T, soil_values::NamedTuple, args...
) where {FT<:AbstractFloat,T<:ModelDimension}
    return initialize(
        FT,
        OSwater,
        Dict{String,Any}(),
        (FT, dimension_value(dim), soil_values.roof.ms, soil_values.ground.ms),
        soil_values,
        args...,
    )
end

"""
    Qinlat{FT<:AbstractFloat, N, MG} <: AbstractLayeredSoilVariables{FT,N}

Lateral soil water flux variables.

# Fields
- `Qin_bare2imp`: Lateral water flux from bare soil to impervious areas [mm/h]
- `Qin_veg2imp`: Lateral water flux from vegetated soil to impervious areas [mm/h]
- `Qin_veg2bare`: Lateral water flux from vegetated soil to bare soil [mm/h]
- `Qin_imp2bare`: Lateral water flux from impervious areas to bare soil [mm/h]
- `Qin_bare2veg`: Lateral water flux from bare soil to vegetated soil [mm/h]
- `Qin_imp2veg`: Lateral water flux from impervious areas to vegetated soil [mm/h]
- `Qin_imp`: Total lateral water flux to impervious areas [mm/h]
- `Qin_bare`: Total lateral water flux to bare soil [mm/h]
- `Qin_veg`: Total lateral water flux to vegetated soil [mm/h]
"""
Base.@kwdef mutable struct Qinlat{FT<:AbstractFloat,N,MG} <:
                           AbstractLayeredSoilVariables{FT,N}
    Qin_bare2imp::Array{MVector{MG,FT},N}
    Qin_veg2imp::Array{MVector{MG,FT},N}
    Qin_veg2bare::Array{MVector{MG,FT},N}
    Qin_imp2bare::Array{MVector{MG,FT},N}
    Qin_bare2veg::Array{MVector{MG,FT},N}
    Qin_imp2veg::Array{MVector{MG,FT},N}
    Qin_imp::Array{MVector{MG,FT},N}
    Qin_bare::Array{MVector{MG,FT},N}
    Qin_veg::Array{MVector{MG,FT},N}
end

function ground_fields(::Type{Qinlat})
    return [
        "Qin_bare2imp",
        "Qin_veg2imp",
        "Qin_veg2bare",
        "Qin_imp2bare",
        "Qin_bare2veg",
        "Qin_imp2veg",
        "Qin_imp",
        "Qin_bare",
        "Qin_veg",
    ]
end

function initialize_qinlat(
    ::Type{FT}, dim::T, soil_values::NamedTuple, args...
) where {FT<:AbstractFloat,T<:ModelDimension}
    return initialize(
        FT,
        Qinlat,
        Dict{String,Any}(),
        (FT, dimension_value(dim), soil_values.ground.ms),
        soil_values,
        args...,
    )
end

"""
    ExWater{FT<:AbstractFloat, N, M} <: AbstractLayeredSoilVariables{FT,N}

Extractable water for plants from soil.

# Fields
- `ExWaterRoofVeg_H`: Extractable water for high vegetation from roof soil [mm m²/m² ground h]
- `ExWaterRoofVeg_L`: Extractable water for low vegetation from roof soil [mm m²/m² ground h]
- `ExWaterGroundImp_H`: Extractable water for high vegetation from impervious ground soil [mm m²/m² ground h]
- `ExWaterGroundImp_L`: Extractable water for low vegetation from impervious ground soil [mm m²/m² ground h]
- `ExWaterGroundBare_H`: Extractable water for high vegetation from bare ground soil [mm m²/m² ground h]
- `ExWaterGroundBare_L`: Extractable water for low vegetation from bare ground soil [mm m²/m² ground h]
- `ExWaterGroundVeg_H`: Extractable water for high vegetation from vegetated ground soil [mm m²/m² ground h]
- `ExWaterGroundVeg_L`: Extractable water for low vegetation from vegetated ground soil [mm m²/m² ground h]
- `ExWaterGroundTot_H`: Extractable water for high vegetation from total ground soil [mm m²/m² ground h]
- `ExWaterGroundTot_L`: Extractable water for low vegetation from total ground soil [mm m²/m² ground h]
"""
# Same names as SoilPotW
Base.@kwdef mutable struct ExWater{FT<:AbstractFloat,N,MR,MG} <:
                           AbstractLayeredSoilVariables{FT,N}
    ExWaterRoofVeg_H::Array{MVector{MR,FT},N}
    ExWaterRoofVeg_L::Array{MVector{MR,FT},N}
    ExWaterGroundImp_H::Array{MVector{MG,FT},N}
    ExWaterGroundImp_L::Array{MVector{MG,FT},N}
    ExWaterGroundBare_H::Array{MVector{MG,FT},N}
    ExWaterGroundBare_L::Array{MVector{MG,FT},N}
    ExWaterGroundVeg_H::Array{MVector{MG,FT},N}
    ExWaterGroundVeg_L::Array{MVector{MG,FT},N}
    ExWaterGroundTot_H::Array{MVector{MG,FT},N}
    ExWaterGroundTot_L::Array{MVector{MG,FT},N}
end

function roof_fields(::Type{ExWater})
    return ["ExWaterRoofVeg_H", "ExWaterRoofVeg_L"]
end

function ground_fields(::Type{ExWater})
    return [
        "ExWaterGroundImp_H",
        "ExWaterGroundImp_L",
        "ExWaterGroundBare_H",
        "ExWaterGroundBare_L",
        "ExWaterGroundVeg_H",
        "ExWaterGroundVeg_L",
        "ExWaterGroundTot_H",
        "ExWaterGroundTot_L",
    ]
end

function initialize_exwater(
    ::Type{FT}, dim::T, soil_values::NamedTuple, args...
) where {FT<:AbstractFloat,T<:ModelDimension}
    return initialize(
        FT,
        ExWater,
        Dict{String,Any}(),
        (FT, dimension_value(dim), soil_values.roof.ms, soil_values.ground.ms),
        soil_values,
        args...,
    )
end

"""
    SoilPotW{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Soil water potential for plants.

# Fields
- `SoilPotWRoofVeg_H`: Soil water potential for plants, high roof vegetation [MPa]
- `SoilPotWRoofVeg_L`: Soil water potential for plants, low roof vegetation [MPa]
- `SoilPotWGroundImp_H`: Soil water potential for plants in impervious ground, high vegetation [MPa]
- `SoilPotWGroundImp_L`: Soil water potential for plants in impervious ground, low vegetation [MPa]
- `SoilPotWGroundBare_H`: Soil water potential for plants in bare ground, high vegetation [MPa]
- `SoilPotWGroundBare_L`: Soil water potential for plants in bare ground, low vegetation [MPa]
- `SoilPotWGroundVeg_H`: Soil water potential for plants in vegetated ground, high vegetation [MPa]
- `SoilPotWGroundVeg_L`: Soil water potential for plants in vegetated ground, low vegetation [MPa]
- `SoilPotWGroundTot_H`: Soil water potential for plants in total ground, high vegetation [MPa]
- `SoilPotWGroundTot_L`: Soil water potential for plants in total ground, low vegetation [MPa]
"""
Base.@kwdef mutable struct SoilPotW{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    SoilPotWRoofVeg_H::Array{FT,N}
    SoilPotWRoofVeg_L::Array{FT,N}
    SoilPotWGroundImp_H::Array{FT,N}
    SoilPotWGroundImp_L::Array{FT,N}
    SoilPotWGroundBare_H::Array{FT,N}
    SoilPotWGroundBare_L::Array{FT,N}
    SoilPotWGroundVeg_H::Array{FT,N}
    SoilPotWGroundVeg_L::Array{FT,N}
    SoilPotWGroundTot_H::Array{FT,N}
    SoilPotWGroundTot_L::Array{FT,N}
end

function initialize_soilpotw(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, SoilPotW, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_soilpotw(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, SoilPotW, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    CiCO2Leaf{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Intercellular CO2 concentration in leaf for different urban surfaces.

# Fields
- `CiCO2LeafRoofVegSun`: Ci_sun_veg sunlit roof leafs [umolCO2/mol]
- `CiCO2LeafRoofVegShd`: Ci_shd_veg shaded roof leafs [umolCO2/mol]
- `CiCO2LeafGroundVegSun`: Ci_sun_veg sunlit ground leafs [umolCO2/mol]
- `CiCO2LeafGroundVegShd`: Ci_shd_veg shaded ground leafs [umolCO2/mol]
- `CiCO2LeafTreeSun`: Ci_sun_veg sunlit tree leafs [umolCO2/mol]
- `CiCO2LeafTreeShd`: Ci_shd_veg shaded tree leafs [umolCO2/mol]
"""
Base.@kwdef mutable struct CiCO2Leaf{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    CiCO2LeafRoofVegSun::Array{FT,N}
    CiCO2LeafRoofVegShd::Array{FT,N}
    CiCO2LeafGroundVegSun::Array{FT,N}
    CiCO2LeafGroundVegShd::Array{FT,N}
    CiCO2LeafTreeSun::Array{FT,N}
    CiCO2LeafTreeShd::Array{FT,N}
end

function initialize_cico2leaf(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, CiCO2Leaf, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_cico2leaf(
    ::Type{FT}, ::TimeSeries, initial_value::FT, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        CiCO2Leaf,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
        initial_value,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{CiCO2Leaf},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    initial_value::FT,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(CiCO2Leaf, data, params, hours)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    field_names = fieldnames(CiCO2Leaf)

    if params[2] == 1
        for var in field_names
            processed[String(var)][1] = initial_value
        end
    end

    return processed
end

"""
    WaterFluxVariables{FT<:AbstractFloat, N} <: Abstract1PModelVariablesSet{FT, N}

Container for all water flux variable components.

# Fields
- `Eflux`: Evapotranspiration fluxes
- `Runoff`: Runoff variables for different urban surfaces
- `Runon`: Runon variables for urban area
- `Leakage`: Leakage variables for different urban surfaces
- `Interception`: Interception variables for different urban surfaces
- `dInt_dt`: Change in interception variables
- `Infiltration`: Infiltration variables for different urban surfaces
- `Vwater`: Water volume in soil for different urban surfaces
- `dVwater_dt`: Change in water volume in soil
- `Owater`: Soil moisture in different soil layers
- `OSwater`: Additional soil moisture variables
- `Qinlat`: Lateral soil water flux variables
- `ExWater`: Extractable water for plants from soil
- `SoilPotW`: Soil water potential for plants
- `CiCO2Leaf`: Intercellular CO2 concentration in leaf for different urban surfaces
"""
Base.@kwdef struct WaterFluxVariables{FT<:AbstractFloat,N,MR,MG} <:
                   Abstract1PModelVariablesSet{FT,N}
    Eflux::Eflux{FT,N}
    Runoff::Runoff{FT,N}
    Runon::Runon{FT,N}
    Leakage::Leakage{FT,N}
    Interception::Interception{FT,N}
    dInt_dt::dInt_dt{FT,N}
    Infiltration::Infiltration{FT,N}
    Vwater::Vwater{FT,N,MR,MG}
    dVwater_dt::dVwater_dt{FT,N,MR,MG}
    Owater::Owater{FT,N,MR,MG}
    OSwater::OSwater{FT,N,MR,MG}
    Qinlat::Qinlat{FT,N,MG}
    ExWater::ExWater{FT,N,MR,MG}
    SoilPotW::SoilPotW{FT,N}
    CiCO2Leaf::CiCO2Leaf{FT,N}
end

function initialize_water_flux_variables(
    ::Type{FT},
    dim::T,
    soil_parameters::SoilParameters{FT},
    vegetation_parameters::VegetationParameters{FT},
    args...,
) where {FT<:AbstractFloat,T<:ModelDimension}
    N = dimension_value(dim)
    return initialize(
        FT,
        WaterFluxVariables,
        Dict{String,Any}(),
        (FT, N, soil_parameters.roof.ms, soil_parameters.ground.ms),
        soil_parameters,
        vegetation_parameters,
        args...,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{WaterFluxVariables},
    data::Dict{String,Any},
    params::Tuple,
    soil_parameters::SoilParameters{FT},
    vegetation_parameters::VegetationParameters{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    roof_soil_params = calculate_soil_values(
        soil_parameters.roof, vegetation_parameters.roof, vegetation_parameters.roof
    )
    ground_soil_params = calculate_soil_values(
        soil_parameters.ground, vegetation_parameters.tree, vegetation_parameters.ground
    )

    soil_values = (; roof=roof_soil_params, ground=ground_soil_params)

    processed["Eflux"] = initialize_eflux(FT, dimensionality_type(params[2]))
    processed["Runoff"] = initialize_runoff(FT, dimensionality_type(params[2]))
    processed["Runon"] = initialize_runon(FT, dimensionality_type(params[2]))
    processed["Leakage"] = initialize_leakage(FT, dimensionality_type(params[2]))
    processed["Interception"] = initialize_interception(FT, dimensionality_type(params[2]))
    processed["dInt_dt"] = initialize_dint_dt(FT, dimensionality_type(params[2]))
    processed["Infiltration"] = initialize_infiltration(FT, dimensionality_type(params[2]))
    processed["Vwater"] = initialize_vwater(FT, dimensionality_type(params[2]), soil_values)
    processed["dVwater_dt"] = initialize_dvwater_dt(
        FT, dimensionality_type(params[2]), soil_values
    )
    processed["Owater"] = initialize_owater(FT, dimensionality_type(params[2]), soil_values)
    processed["OSwater"] = initialize_oswater(
        FT, dimensionality_type(params[2]), soil_values
    )
    processed["Qinlat"] = initialize_qinlat(FT, dimensionality_type(params[2]), soil_values)
    processed["ExWater"] = initialize_exwater(
        FT, dimensionality_type(params[2]), soil_values
    )
    processed["SoilPotW"] = initialize_soilpotw(FT, dimensionality_type(params[2]))
    processed["CiCO2Leaf"] = initialize_cico2leaf(FT, dimensionality_type(params[2]))

    return processed
end

function calculate_soil_values(
    soil::VegetatedSoilParameters{FT},
    tree::HeightDependentVegetationParameters{FT},
    ground::HeightDependentVegetationParameters{FT},
) where {FT<:AbstractFloat}
    _, _, _, Osat, Ohy, _, _, _, _, _, O33 = soil_parameters_total(
        soil.Pcla,
        soil.Psan,
        soil.Porg,
        soil.Kfc,
        soil.Phy,
        soil.SPAR,
        soil.Kbot,
        tree.CASE_ROOT,
        ground.CASE_ROOT,
        tree.ZR95,
        ground.ZR95,
        tree.ZR50,
        ground.ZR50,
        tree.ZRmax,
        ground.ZRmax,
        soil.Zs,
    )

    unique_Osat = unique(Osat)
    @assert length(unique_Osat) == 1 "Osat should be unique after preprocessing"
    unique_Ohy = unique(Ohy)
    @assert length(unique_Ohy) == 1 "Ohy should be unique after preprocessing"
    unique_O33 = unique(O33)
    @assert length(unique_O33) == 1 "O33 should be unique after preprocessing"

    return (;
        Osat=unique_Osat[], Ohy=unique_Ohy[], O33=unique_O33[], ms=soil.ms, dz=diff(soil.Zs)
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{WaterFluxVariables},
    data::Dict{String,Any},
    params::Tuple,
    soil_parameters::SoilParameters{FT},
    vegetation_parameters::VegetationParameters{FT},
    initial_value::FT,
    hours::Int,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    roof_soil_params = calculate_soil_values(
        soil_parameters.roof, vegetation_parameters.roof, vegetation_parameters.roof
    )
    ground_soil_params = calculate_soil_values(
        soil_parameters.ground, vegetation_parameters.tree, vegetation_parameters.ground
    )

    soil_values = (; roof=roof_soil_params, ground=ground_soil_params)

    processed["Eflux"] = initialize_eflux(FT, dimensionality_type(params[2]), hours)
    processed["Runoff"] = initialize_runoff(FT, dimensionality_type(params[2]), hours)
    processed["Runon"] = initialize_runon(FT, dimensionality_type(params[2]), hours)
    processed["Leakage"] = initialize_leakage(FT, dimensionality_type(params[2]), hours)
    processed["Interception"] = initialize_interception(
        FT, dimensionality_type(params[2]), hours
    )
    processed["dInt_dt"] = initialize_dint_dt(FT, dimensionality_type(params[2]), hours)
    processed["Infiltration"] = initialize_infiltration(
        FT, dimensionality_type(params[2]), hours
    )
    processed["Vwater"] = initialize_vwater(
        FT, dimensionality_type(params[2]), soil_values, hours
    )
    processed["dVwater_dt"] = initialize_dvwater_dt(
        FT, dimensionality_type(params[2]), soil_values, hours
    )
    processed["Owater"] = initialize_owater(
        FT, dimensionality_type(params[2]), soil_values, hours
    )
    processed["OSwater"] = initialize_oswater(
        FT, dimensionality_type(params[2]), soil_values, hours
    )
    processed["Qinlat"] = initialize_qinlat(
        FT, dimensionality_type(params[2]), soil_values, hours
    )
    processed["ExWater"] = initialize_exwater(
        FT, dimensionality_type(params[2]), soil_values, hours
    )
    processed["SoilPotW"] = initialize_soilpotw(FT, dimensionality_type(params[2]), hours)
    processed["CiCO2Leaf"] = initialize_cico2leaf(
        FT, dimensionality_type(params[2]), initial_value, hours
    )

    return processed
end
