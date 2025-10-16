# Subtype for the set of subsets of water flux model variables
abstract type AbstractWaterFluxVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT} end

# Subtype for each subset of water flux model variables
abstract type AbstractWaterFluxVariablesSubset{FT<:AbstractFloat,N} <:
              AbstractModelVariables{FT} end

abstract type AbstractEflux{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractRunoff{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractRunon{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractLeakage{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractInterception{FT<:AbstractFloat,N} <:
              AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractdInt_dt{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractInfiltration{FT<:AbstractFloat,N} <:
              AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractVwater{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractdVwater_dt{FT<:AbstractFloat,N} <:
              AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractOwater{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractOSwater{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractQinlat{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractSoilPotW{FT<:AbstractFloat,N} <:
              AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractExWater{FT<:AbstractFloat,N} <: AbstractWaterFluxVariablesSubset{FT,N} end
abstract type AbstractCiCO2Leaf{FT<:AbstractFloat,N} <:
              AbstractWaterFluxVariablesSubset{FT,N} end

"""
    Eflux{FT<:AbstractFloat, N} <: AbstractEflux{FT,N}

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
Base.@kwdef struct Eflux{FT<:AbstractFloat,N} <: AbstractEflux{FT,N}
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

function initialize_eflux(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, Eflux, Dict{String,Any}(), (FT, N), hours)
end

"""
    Runoff{FT<:AbstractFloat, N} <: AbstractRunoff{FT,N}

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
Base.@kwdef struct Runoff{FT<:AbstractFloat,N} <: AbstractRunoff{FT,N}
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

function initialize_runoff(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, Runoff, Dict{String,Any}(), (FT, N), hours)
end

"""
    Runon{FT<:AbstractFloat, N} <: AbstractRunon{FT,N}

Runon variables for urban area.

# Fields
- `RunonRoofTot`: Total roof runon to the next time step [mm/time step per horizontal roof area]
- `RunoffRoofTot`: Total roof runoff that is removed from the system [mm/time step per horizontal roof area]
- `RunonGroundTot`: Total runon in canyon to the next time step [mm/time step per horizontal ground area]
- `RunoffGroundTot`: Total runoff in canyon that is removed from the system [mm/time step per horizontal ground area]
- `RunonUrban`: Total urban runon to the next time step [mm/time step per horizontal urban area]
- `RunoffUrban`: Total urban runoff that is removed from the system [mm/time step per horizontal urban area]
"""
Base.@kwdef struct Runon{FT<:AbstractFloat,N} <: AbstractRunon{FT,N}
    RunonRoofTot::Array{FT,N}
    RunoffRoofTot::Array{FT,N}
    RunonGroundTot::Array{FT,N}
    RunoffGroundTot::Array{FT,N}
    RunonUrban::Array{FT,N}
    RunoffUrban::Array{FT,N}
end

function initialize_runon(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, Runon, Dict{String,Any}(), (FT, N), hours)
end

"""
    Leakage{FT<:AbstractFloat, N} <: AbstractLeakage{FT,N}

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
Base.@kwdef struct Leakage{FT<:AbstractFloat,N} <: AbstractLeakage{FT,N}
    LkRoofImp::Array{FT,N}
    LkRoofVeg::Array{FT,N}
    LkRoof::Array{FT,N}
    LkGroundImp::Array{FT,N}
    LkGroundBare::Array{FT,N}
    LkGroundVeg::Array{FT,N}
    LkGround::Array{FT,N}
    LkUrban::Array{FT,N}
end

function initialize_leakage(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, Leakage, Dict{String,Any}(), (FT, N), hours)
end

"""
    Interception{FT<:AbstractFloat, N} <: AbstractInterception{FT,N}

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
Base.@kwdef struct Interception{FT<:AbstractFloat,N} <: AbstractInterception{FT,N}
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

function initialize_interception(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, Interception, Dict{String,Any}(), (FT, N), hours)
end

"""
    dInt_dt{FT<:AbstractFloat, N} <: AbstractdInt_dt{FT,N}

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
Base.@kwdef struct dInt_dt{FT<:AbstractFloat,N} <: AbstractdInt_dt{FT,N}
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

function initialize_dint_dt(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, dInt_dt, Dict{String,Any}(), (FT, N), hours)
end

"""
    Infiltration{FT<:AbstractFloat, N} <: AbstractInfiltration{FT,N}

Infiltration variables for different urban surfaces.

# Fields
- `fRoofVeg`: Infiltration in first soil layer of vegetated roof [mm/h per horizontal roof area]
- `fGroundBare`: Infiltration in first soil layer of bare ground [mm/h per horizontal bare ground area]
- `fGroundVeg`: Infiltration in first soil layer of vegetated ground [mm/h per horizontal vegetated ground area]
- `fGroundImp`: Infiltration in impervious ground (usually zero) [mm/h per horizontal impervious ground area]
"""
Base.@kwdef struct Infiltration{FT<:AbstractFloat,N} <: AbstractInfiltration{FT,N}
    fRoofVeg::Array{FT,N}
    fGroundBare::Array{FT,N}
    fGroundVeg::Array{FT,N}
    fGroundImp::Array{FT,N}
end

function initialize_infiltration(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, Infiltration, Dict{String,Any}(), (FT, N), hours)
end

"""
    Vwater{FT<:AbstractFloat, N} <: AbstractVwater{FT,N}

Water volume in soil for different urban surfaces.

# Fields
- `VRoofSoilVeg`: Water volume in the different soil layers of roof [mm per horizontal roof area]
- `VGroundSoilImp`: Water volume in the different soil layers of ground under impervious [mm per horizontal impervious ground area]
- `VGroundSoilBare`: Water volume in the different soil layers of ground under bare [mm per horizontal bare ground area]
- `VGroundSoilVeg`: Water volume in the different soil layers of ground under vegetated [mm per horizontal vegetated ground area]
- `VGroundSoilTot`: Water volume in the different soil layers of ground total [mm per horizontal ground area]
"""
Base.@kwdef struct Vwater{FT<:AbstractFloat,N} <: AbstractVwater{FT,N}
    VRoofSoilVeg::Array{FT,N}
    VGroundSoilImp::Array{FT,N}
    VGroundSoilBare::Array{FT,N}
    VGroundSoilVeg::Array{FT,N}
    VGroundSoilTot::Array{FT,N}
end

function initialize_vwater(
    ::Type{FT}, N::Int, soil_parameters::SoilParameters{FT}, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, Vwater, Dict{String,Any}(), (FT, N), hours, soil_parameters)
end

function get_dimensions(
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {T<:AbstractVwater,FT<:AbstractFloat}
    if params[2] ∉ [1, 2]
        throw(ArgumentError("Only N=1 and N=2 are currently supported"))
    end

    roof_layers = soil_parameters.roof.ms
    ground_layers = soil_parameters.ground.ms

    if params[2] == 1
        return Dict{String,Tuple}(
            "VRoofSoilVeg" => (roof_layers,),
            "VGroundSoilImp" => (ground_layers,),
            "VGroundSoilBare" => (ground_layers,),
            "VGroundSoilVeg" => (ground_layers,),
            "VGroundSoilTot" => (ground_layers,),
        )
    else
        return Dict{String,Tuple}(
            "VRoofSoilVeg" => (hours, roof_layers),
            "VGroundSoilImp" => (hours, ground_layers),
            "VGroundSoilBare" => (hours, ground_layers),
            "VGroundSoilVeg" => (hours, ground_layers),
            "VGroundSoilTot" => (hours, ground_layers),
        )
    end
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{Vwater},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(Vwater, data, params, hours, soil_parameters)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    # Initialize temperature and humidity fields
    ground_fields = [
        "VGroundSoilImp", "VGroundSoilBare", "VGroundSoilVeg", "VGroundSoilTot"
    ]
    ground_init = soil_parameters.ground.O33 * soil_parameters.ground.dz

    roof_fields = ["VRoofSoilVeg"]
    roof_init = soil_parameters.roof.O33 * soil_parameters.roof.dz

    if params[2] == 2
        for var in ground_fields
            processed[var][1, :] = ground_init
        end
        for var in roof_fields
            processed[var][1, :] = roof_init
        end
    end

    return processed
end

"""
    dVwater_dt{FT<:AbstractFloat, N} <: AbstractdVwater_dt{FT,N}

Change in water volume in soil for different urban surfaces.

# Fields
- `dVRoofSoilVeg_dt`: Change in water volume in the different soil layers of roof [mm per horizontal roof area]
- `dVGroundSoilImp_dt`: Change in water volume in the different soil layers of ground under impervious [mm per horizontal impervious ground area]
- `dVGroundSoilBare_dt`: Change in water volume in the different soil layers of ground under bare [mm per horizontal bare ground area]
- `dVGroundSoilVeg_dt`: Change in water volume in the different soil layers of ground under vegetated [mm per horizontal ground area]
- `dVGroundSoilTot_dt`: Change in water volume in the different soil layers of ground total [mm per horizontal ground area]
"""
# Same names as Vwater
Base.@kwdef struct dVwater_dt{FT<:AbstractFloat,N} <: AbstractdVwater_dt{FT,N}
    dVRoofSoilVeg_dt::Array{FT,N}
    dVGroundSoilImp_dt::Array{FT,N}
    dVGroundSoilBare_dt::Array{FT,N}
    dVGroundSoilVeg_dt::Array{FT,N}
    dVGroundSoilTot_dt::Array{FT,N}
end

function initialize_dvwater_dt(
    ::Type{FT}, N::Int, soil_parameters::SoilParameters{FT}, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, dVwater_dt, Dict{String,Any}(), (FT, N), hours, soil_parameters)
end

function get_dimensions(
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {T<:AbstractdVwater_dt,FT<:AbstractFloat}
    if params[2] ∉ [1, 2]
        throw(ArgumentError("Only N=1 and N=2 are currently supported"))
    end

    roof_layers = soil_parameters.roof.ms
    ground_layers = soil_parameters.ground.ms

    if params[2] == 1
        return Dict{String,Tuple}(
            "dVRoofSoilVeg_dt" => (roof_layers,),
            "dVGroundSoilImp_dt" => (ground_layers,),
            "dVGroundSoilBare_dt" => (ground_layers,),
            "dVGroundSoilVeg_dt" => (ground_layers,),
            "dVGroundSoilTot_dt" => (ground_layers,),
        )
    else
        return Dict{String,Tuple}(
            "dVRoofSoilVeg_dt" => (hours, roof_layers),
            "dVGroundSoilImp_dt" => (hours, ground_layers),
            "dVGroundSoilBare_dt" => (hours, ground_layers),
            "dVGroundSoilVeg_dt" => (hours, ground_layers),
            "dVGroundSoilTot_dt" => (hours, ground_layers),
        )
    end
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {FT<:AbstractFloat,T<:AbstractWaterFluxVariablesSubset}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(T, data, params, hours, soil_parameters)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    return processed
end

"""
    Owater{FT<:AbstractFloat, N} <: AbstractOwater{FT,N}

Soil moisture in different soil layers for urban surfaces.

# Fields
- `OwRoofSoilVeg`: Soil moisture in the different soil layers of roof [-]
- `OwGroundSoilImp`: Soil moisture in the different soil layers of ground under impervious [-]
- `OwGroundSoilBare`: Soil moisture in the different soil layers of ground under bare [-]
- `OwGroundSoilVeg`: Soil moisture in the different soil layers of ground under vegetated [-]
- `OwGroundSoilTot`: Soil moisture in the different soil layers of ground total [-]
"""
# Same names as Vwater
Base.@kwdef struct Owater{FT<:AbstractFloat,N} <: AbstractOwater{FT,N}
    OwRoofSoilVeg::Array{FT,N}
    OwGroundSoilImp::Array{FT,N}
    OwGroundSoilBare::Array{FT,N}
    OwGroundSoilVeg::Array{FT,N}
    OwGroundSoilTot::Array{FT,N}
end

function initialize_owater(
    ::Type{FT}, N::Int, soil_parameters::SoilParameters{FT}, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, Owater, Dict{String,Any}(), (FT, N), hours, soil_parameters)
end

function get_dimensions(
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {T<:AbstractOwater,FT<:AbstractFloat}
    if params[2] ∉ [1, 2]
        throw(ArgumentError("Only N=1 and N=2 are currently supported"))
    end

    roof_layers = soil_parameters.roof.ms
    ground_layers = soil_parameters.ground.ms

    if params[2] == 1
        return Dict{String,Tuple}(
            "OwRoofSoilVeg" => (roof_layers,),
            "OwGroundSoilImp" => (ground_layers,),
            "OwGroundSoilBare" => (ground_layers,),
            "OwGroundSoilVeg" => (ground_layers,),
            "OwGroundSoilTot" => (ground_layers,),
        )
    else
        return Dict{String,Tuple}(
            "OwRoofSoilVeg" => (hours, roof_layers),
            "OwGroundSoilImp" => (hours, ground_layers),
            "OwGroundSoilBare" => (hours, ground_layers),
            "OwGroundSoilVeg" => (hours, ground_layers),
            "OwGroundSoilTot" => (hours, ground_layers),
        )
    end
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{Owater},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(Owater, data, params, hours, soil_parameters)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    # Initialize temperature and humidity fields
    ground_fields = [
        "OwGroundSoilImp", "OwGroundSoilBare", "OwGroundSoilVeg", "OwGroundSoilTot"
    ]
    roof_fields = ["OwRoofSoilVeg"]

    if params[2] == 2
        for var in ground_fields
            processed[var][1, :] .= soil_parameters.ground.O33
        end
        for var in roof_fields
            processed[var][1, :] .= soil_parameters.roof.O33
        end
    end

    return processed
end

"""
    OSwater{FT<:AbstractFloat, N} <: AbstractOSwater{FT,N}

Additional soil moisture variables for urban surfaces.

# Fields
- `OSwRoofSoilVeg`: Additional soil moisture values for roof soil layers [-]
- `OSwGroundSoilImp`: Additional soil moisture values for ground soil layers under impervious [-]
- `OSwGroundSoilBare`: Additional soil moisture values for ground soil layers under bare [-]
- `OSwGroundSoilVeg`: Additional soil moisture values for ground soil layers under vegetated [-]
- `OSwGroundSoilTot`: Additional soil moisture values for ground soil layers total [-]
"""
# Same names as Vwater
Base.@kwdef struct OSwater{FT<:AbstractFloat,N} <: AbstractOSwater{FT,N}
    OSwRoofSoilVeg::Array{FT,N}
    OSwGroundSoilImp::Array{FT,N}
    OSwGroundSoilBare::Array{FT,N}
    OSwGroundSoilVeg::Array{FT,N}
    OSwGroundSoilTot::Array{FT,N}
end

function initialize_oswater(
    ::Type{FT}, N::Int, soil_parameters::SoilParameters{FT}, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, OSwater, Dict{String,Any}(), (FT, N), hours, soil_parameters)
end

function get_dimensions(
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {T<:AbstractOSwater,FT<:AbstractFloat}
    if params[2] ∉ [1, 2]
        throw(ArgumentError("Only N=1 and N=2 are currently supported"))
    end

    roof_layers = soil_parameters.roof.ms
    ground_layers = soil_parameters.ground.ms

    if params[2] == 1
        return Dict{String,Tuple}(
            "OSwRoofSoilVeg" => (roof_layers,),
            "OSwGroundSoilImp" => (ground_layers,),
            "OSwGroundSoilBare" => (ground_layers,),
            "OSwGroundSoilVeg" => (ground_layers,),
            "OSwGroundSoilTot" => (ground_layers,),
        )
    else
        return Dict{String,Tuple}(
            "OSwRoofSoilVeg" => (hours, roof_layers),
            "OSwGroundSoilImp" => (hours, ground_layers),
            "OSwGroundSoilBare" => (hours, ground_layers),
            "OSwGroundSoilVeg" => (hours, ground_layers),
            "OSwGroundSoilTot" => (hours, ground_layers),
        )
    end
end

"""
    Qinlat{FT<:AbstractFloat, N, M} <: AbstractQinlat{FT,N}

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
Base.@kwdef struct Qinlat{FT<:AbstractFloat,N} <: AbstractQinlat{FT,N}
    Qin_bare2imp::Array{FT,N}
    Qin_veg2imp::Array{FT,N}
    Qin_veg2bare::Array{FT,N}
    Qin_imp2bare::Array{FT,N}
    Qin_bare2veg::Array{FT,N}
    Qin_imp2veg::Array{FT,N}
    Qin_imp::Array{FT,N}
    Qin_bare::Array{FT,N}
    Qin_veg::Array{FT,N}
end

function initialize_qinlat(
    ::Type{FT}, N::Int, soil_parameters::SoilParameters{FT}, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(
        FT, Qinlat, Dict{String,Any}(), (FT, N), hours, soil_parameters.ground.ms
    )
end

function get_dimensions(
    ::Type{T}, data::Dict{String,Any}, params::Tuple{DataType,Signed}, hours::Int, ms::Int
) where {T<:AbstractQinlat}
    if params[2] ∉ [1, 2]
        throw(ArgumentError("Only N=1 and N=2 are currently supported"))
    end

    # Create a dictionary with all field names and their dimensions
    dimensions = Dict{String,Tuple}()

    # Get all field names from the struct
    field_names = fieldnames(T)

    if params[2] == 1
        # For scalar case (N=1), all fields are vectors with length ms
        for field in field_names
            dimensions[String(field)] = (ms,)
        end
    else
        # For time series case (N=2), all fields are matrices with size (hours, ms)
        for field in field_names
            dimensions[String(field)] = (hours, ms)
        end
    end

    return dimensions
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{T}, data::Dict{String,Any}, params::Tuple, hours::Int, ms::Int
) where {FT<:AbstractFloat,T<:AbstractQinlat}
    processed = Dict{String,Any}()
    dimensions = get_dimensions(T, data, params, hours, ms)

    for (var, dims) in dimensions
        processed[var] = zeros(FT, dims)
    end

    return processed
end

"""
    ExWater{FT<:AbstractFloat, N, M} <: AbstractExWater{FT,N}

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
Base.@kwdef struct ExWater{FT<:AbstractFloat,N} <: AbstractExWater{FT,N}
    ExWaterRoofVeg_H::Array{FT,N}
    ExWaterRoofVeg_L::Array{FT,N}
    ExWaterGroundImp_H::Array{FT,N}
    ExWaterGroundImp_L::Array{FT,N}
    ExWaterGroundBare_H::Array{FT,N}
    ExWaterGroundBare_L::Array{FT,N}
    ExWaterGroundVeg_H::Array{FT,N}
    ExWaterGroundVeg_L::Array{FT,N}
    ExWaterGroundTot_H::Array{FT,N}
    ExWaterGroundTot_L::Array{FT,N}
end

function initialize_exwater(
    ::Type{FT}, N::Int, soil_parameters::SoilParameters{FT}, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, ExWater, Dict{String,Any}(), (FT, N), hours, soil_parameters)
end

function get_dimensions(
    ::Type{T},
    data::Dict{String,Any},
    params::Tuple{DataType,Signed},
    hours::Int,
    soil_parameters::SoilParameters{FT},
) where {T<:AbstractExWater,FT<:AbstractFloat}
    if params[2] ∉ [1, 2]
        throw(ArgumentError("Only N=1 and N=2 are currently supported"))
    end

    ms_ground = soil_parameters.ground.ms
    ms_roof = soil_parameters.roof.ms

    ground_names = [
        "ExWaterGroundImp_H",
        "ExWaterGroundImp_L",
        "ExWaterGroundBare_H",
        "ExWaterGroundBare_L",
        "ExWaterGroundVeg_H",
        "ExWaterGroundVeg_L",
        "ExWaterGroundTot_H",
        "ExWaterGroundTot_L",
    ]
    roof_names = ["ExWaterRoofVeg_H", "ExWaterRoofVeg_L"]

    # Create a dictionary with all field names and their dimensions
    dimensions = Dict{String,Tuple}()

    if params[2] == 1
        # For scalar case (N=1), all fields are vectors with length ms
        for field in roof_names
            dimensions[field] = (ms_roof,)
        end
        for field in ground_names
            dimensions[field] = (ms_ground,)
        end
    else
        # For time series case (N=2), all fields are matrices with size (hours, ms)
        for field in roof_names
            dimensions[field] = (hours, ms_roof)
        end
        for field in ground_names
            dimensions[field] = (hours, ms_ground)
        end
    end

    return dimensions
end

"""
    SoilPotW{FT<:AbstractFloat, N} <: AbstractSoilPotW{FT,N}

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
Base.@kwdef struct SoilPotW{FT<:AbstractFloat,N} <: AbstractSoilPotW{FT,N}
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
function initialize_soilpotw(::Type{FT}, N::Int, hours::Int=1) where {FT<:AbstractFloat}
    return initialize(FT, SoilPotW, Dict{String,Any}(), (FT, N), hours)
end

"""
    CiCO2Leaf{FT<:AbstractFloat, N} <: AbstractCiCO2Leaf{FT,N}

Intercellular CO2 concentration in leaf for different urban surfaces.

# Fields
- `CiCO2LeafRoofVegSun`: Ci_sun_veg sunlit roof leafs [umolCO2/mol]
- `CiCO2LeafRoofVegShd`: Ci_shd_veg shaded roof leafs [umolCO2/mol]
- `CiCO2LeafGroundVegSun`: Ci_sun_veg sunlit ground leafs [umolCO2/mol]
- `CiCO2LeafGroundVegShd`: Ci_shd_veg shaded ground leafs [umolCO2/mol]
- `CiCO2LeafTreeSun`: Ci_sun_veg sunlit tree leafs [umolCO2/mol]
- `CiCO2LeafTreeShd`: Ci_shd_veg shaded tree leafs [umolCO2/mol]
"""
Base.@kwdef struct CiCO2Leaf{FT<:AbstractFloat,N} <: AbstractCiCO2Leaf{FT,N}
    CiCO2LeafRoofVegSun::Array{FT,N}
    CiCO2LeafRoofVegShd::Array{FT,N}
    CiCO2LeafGroundVegSun::Array{FT,N}
    CiCO2LeafGroundVegShd::Array{FT,N}
    CiCO2LeafTreeSun::Array{FT,N}
    CiCO2LeafTreeShd::Array{FT,N}
end

function initialize_cico2leaf(
    ::Type{FT}, N::Int, initial_value::FT=400.0, hours::Int=1
) where {FT<:AbstractFloat}
    return initialize(FT, CiCO2Leaf, Dict{String,Any}(), (FT, N), hours, initial_value)
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
    WaterFluxVariables{FT<:AbstractFloat, N, MS} <: AbstractWaterFluxVariables{FT}

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
Base.@kwdef struct WaterFluxVariables{FT<:AbstractFloat,N,Np1} <:
                   AbstractWaterFluxVariables{FT}
    Eflux::Eflux{FT,N}
    Runoff::Runoff{FT,N}
    Runon::Runon{FT,N}
    Leakage::Leakage{FT,N}
    Interception::Interception{FT,N}
    dInt_dt::dInt_dt{FT,N}
    Infiltration::Infiltration{FT,N}
    Vwater::Vwater{FT,Np1}
    dVwater_dt::dVwater_dt{FT,Np1}
    Owater::Owater{FT,Np1}
    OSwater::OSwater{FT,Np1}
    Qinlat::Qinlat{FT,Np1}
    ExWater::ExWater{FT,Np1}
    SoilPotW::SoilPotW{FT,N}
    CiCO2Leaf::CiCO2Leaf{FT,N}
end

# Base getproperty methods for scalar access
function Base.getproperty(
    obj::T, field::Symbol
) where {FT<:AbstractFloat,T<:AbstractWaterFluxVariablesSubset{FT,0}}
    return getfield(obj, field)[]
end

function initialize_water_flux_variables(
    ::Type{FT},
    N::Int,
    soil_parameters::SoilParameters{FT},
    initial_value::FT=400.0,
    hours::Int=1,
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        WaterFluxVariables,
        Dict{String,Any}(),
        (FT, N, N+1),
        hours,
        soil_parameters,
        initial_value,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{WaterFluxVariables},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
    soil_parameters::SoilParameters{FT},
    initial_value::FT,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["Eflux"] = initialize_eflux(FT, params[2], hours)
    processed["Runoff"] = initialize_runoff(FT, params[2], hours)
    processed["Runon"] = initialize_runon(FT, params[2], hours)
    processed["Leakage"] = initialize_leakage(FT, params[2], hours)
    processed["Interception"] = initialize_interception(FT, params[2], hours)
    processed["dInt_dt"] = initialize_dint_dt(FT, params[2], hours)
    processed["Infiltration"] = initialize_infiltration(FT, params[2], hours)
    processed["Vwater"] = initialize_vwater(FT, params[2] + 1, soil_parameters, hours)
    processed["dVwater_dt"] = initialize_dvwater_dt(
        FT, params[2] + 1, soil_parameters, hours
    )
    processed["Owater"] = initialize_owater(FT, params[2] + 1, soil_parameters, hours)
    processed["OSwater"] = initialize_oswater(FT, params[2] + 1, soil_parameters, hours)
    processed["Qinlat"] = initialize_qinlat(FT, params[2] + 1, soil_parameters, hours)
    processed["ExWater"] = initialize_exwater(FT, params[2] + 1, soil_parameters, hours)
    processed["SoilPotW"] = initialize_soilpotw(FT, params[2], hours)
    processed["CiCO2Leaf"] = initialize_cico2leaf(FT, params[2], initial_value, hours)

    return processed
end
