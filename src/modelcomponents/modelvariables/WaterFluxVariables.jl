"""
    Eflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct Eflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    EfluxRoofImp::FT
    EfluxRoofVegInt::FT
    EfluxRoofVegPond::FT
    EfluxRoofVegSoil::FT
    TEfluxRoofVeg::FT
    EfluxRoofVeg::FT
    EfluxRoof::FT
    EfluxGroundImp::FT
    EfluxGroundBarePond::FT
    EfluxGroundBareSoil::FT
    EfluxGroundBare::FT
    EfluxGroundVegInt::FT
    EfluxGroundVegPond::FT
    EfluxGroundVegSoil::FT
    TEfluxGroundVeg::FT
    EfluxGroundVeg::FT
    EfluxGround::FT
    EfluxTreeInt::FT
    TEfluxTree::FT
    EfluxTree::FT
    EfluxWallSun::FT
    EfluxWallShade::FT
    EfluxCanyon::FT
    EfluxUrban::FT
end

function Eflux(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Eflux, Dict{String,Any}())
end

function ground_fields(::Type{Eflux})
    return (
        :EfluxGroundImp,
        :EfluxGroundBarePond,
        :EfluxGroundBareSoil,
        :EfluxGroundBare,
        :EfluxGroundVegInt,
        :EfluxGroundVegPond,
        :EfluxGroundVegSoil,
        :TEfluxGroundVeg,
        :EfluxGroundVeg,
        :EfluxGround,
        :EfluxTreeInt,
        :TEfluxTree,
        :EfluxTree,
        :EfluxWallSun,
        :EfluxWallShade,
        :EfluxCanyon,
    )
end

function roof_fields(::Type{Eflux})
    return (
        :EfluxRoofImp,
        :EfluxRoofVegInt,
        :EfluxRoofVegPond,
        :EfluxRoofVegSoil,
        :TEfluxRoofVeg,
        :EfluxRoofVeg,
        :EfluxRoof,
    )
end

"""
    Runoff{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct Runoff{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    QRoofImp::FT
    QRoofVegDrip::FT
    QRoofVegPond::FT
    QRoofVegSoil::FT
    QGroundImp::FT
    QGroundBarePond::FT
    QGroundBareSoil::FT
    QTree::FT
    QGroundVegDrip::FT
    QGroundVegPond::FT
    QGroundVegSoil::FT
end

function Runoff(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Runoff, Dict{String,Any}())
end

function roof_fields(::Type{Runoff})
    return (:QRoofImp, :QRoofVegDrip, :QRoofVegPond, :QRoofVegSoil)
end

function ground_fields(::Type{Runoff})
    return (
        :QGroundImp,
        :QGroundBarePond,
        :QGroundBareSoil,
        :QTree,
        :QGroundVegDrip,
        :QGroundVegPond,
        :QGroundVegSoil,
    )
end

"""
    Runon{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Runon variables for urban area.

# Fields
- `RunonRoofTot`: Total roof runon to the next time step [mm/time step per horizontal roof area]
- `RunoffRoofTot`: Total roof runoff that is removed from the system [mm/time step per horizontal roof area]
- `RunonGroundTot`: Total runon in canyon to the next time step [mm/time step per horizontal ground area]
- `RunoffGroundTot`: Total runoff in canyon that is removed from the system [mm/time step per horizontal ground area]
- `RunonUrban`: Total urban runon to the next time step [mm/time step per horizontal urban area]
- `RunoffUrban`: Total urban runoff that is removed from the system [mm/time step per horizontal urban area]
"""
Base.@kwdef mutable struct Runon{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    RunonRoofTot::FT
    RunoffRoofTot::FT
    RunonGroundTot::FT
    RunoffGroundTot::FT
    RunonUrban::FT
    RunoffUrban::FT
end

function Runon(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Runon, Dict{String,Any}())
end

function Runon(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return Runon{FT}(;
        RunonRoofTot=data["RunonRoofTot"],
        RunoffRoofTot=data["RunoffRoofTot"],
        RunonGroundTot=data["RunonGroundTot"],
        RunoffGroundTot=data["RunoffGroundTot"],
        RunonUrban=data["RunonUrban"],
        RunoffUrban=data["RunoffUrban"],
    )
end

function update!(dest::Runon{FT}, src::Runon{FT}) where {FT<:AbstractFloat}
    dest.RunonRoofTot = src.RunonRoofTot
    dest.RunoffRoofTot = src.RunoffRoofTot
    dest.RunonGroundTot = src.RunonGroundTot
    dest.RunoffGroundTot = src.RunoffGroundTot
    dest.RunonUrban = src.RunonUrban
    dest.RunoffUrban = src.RunoffUrban

    return nothing
end

function roof_fields(::Type{Runon})
    return (:RunonRoofTot, :RunoffRoofTot)
end

function ground_fields(::Type{Runon})
    return (:RunonGroundTot, :RunoffGroundTot)
end

"""
    Leakage{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct Leakage{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    LkRoofImp::FT
    LkRoofVeg::FT
    LkRoof::FT
    LkGroundImp::FT
    LkGroundBare::FT
    LkGroundVeg::FT
    LkGround::FT
    LkUrban::FT
end

function Leakage(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Leakage, Dict{String,Any}())
end

function roof_fields(::Type{Leakage})
    return (:LkRoofImp, :LkRoofVeg, :LkRoof)
end

function ground_fields(::Type{Leakage})
    return (:LkGroundImp, :LkGroundBare, :LkGroundVeg, :LkGround)
end

"""
    Interception{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct Interception{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    IntRoofImp::FT
    IntRoofVegPlant::FT
    IntRoofVegGround::FT
    IntRooftot::FT
    IntGroundImp::FT
    IntGroundBare::FT
    IntGroundVegPlant::FT
    IntGroundVegGround::FT
    IntTree::FT
end

function Interception(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Interception, Dict{String,Any}())
end

function Interception(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return Interception{FT}(;
        IntRoofImp=data["IntRoofImp"],
        IntRoofVegPlant=data["IntRoofVegPlant"],
        IntRoofVegGround=data["IntRoofVegGround"],
        IntRooftot=data["IntRooftot"],
        IntGroundImp=data["IntGroundImp"],
        IntGroundBare=data["IntGroundBare"],
        IntGroundVegPlant=data["IntGroundVegPlant"],
        IntGroundVegGround=data["IntGroundVegGround"],
        IntTree=data["IntTree"],
    )
end

function update!(dest::Interception{FT}, src::Interception{FT}) where {FT<:AbstractFloat}
    dest.IntRoofImp = src.IntRoofImp
    dest.IntRoofVegPlant = src.IntRoofVegPlant
    dest.IntRoofVegGround = src.IntRoofVegGround
    dest.IntRooftot = src.IntRooftot
    dest.IntGroundImp = src.IntGroundImp
    dest.IntGroundBare = src.IntGroundBare
    dest.IntGroundVegPlant = src.IntGroundVegPlant
    dest.IntGroundVegGround = src.IntGroundVegGround
    dest.IntTree = src.IntTree
    return nothing
end

function roof_fields(::Type{Interception})
    return (:IntRoofImp, :IntRoofVegPlant, :IntRoofVegGround, :IntRooftot)
end

function ground_fields(::Type{Interception})
    return (
        :IntGroundImp, :IntGroundBare, :IntGroundVegPlant, :IntGroundVegGround, :IntTree
    )
end

"""
    dInt_dt{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct dInt_dt{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    dInt_dtRoofImp::FT
    dInt_dtRoofVegPlant::FT
    dInt_dtRoofVegGround::FT
    dInt_dtRooftot::FT
    dInt_dtGroundImp::FT
    dInt_dtGroundBare::FT
    dInt_dtGroundVegPlant::FT
    dInt_dtGroundVegGround::FT
    dInt_dtTree::FT
end

function dInt_dt(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, dInt_dt, Dict{String,Any}())
end

function roof_fields(::Type{dInt_dt})
    return (:dInt_dtRoofImp, :dInt_dtRoofVegPlant, :dInt_dtRoofVegGround, :dInt_dtRooftot)
end

function ground_fields(::Type{dInt_dt})
    return (
        :dInt_dtGroundImp,
        :dInt_dtGroundBare,
        :dInt_dtGroundVegPlant,
        :dInt_dtGroundVegGround,
        :dInt_dtTree,
    )
end

"""
    Infiltration{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Infiltration variables for different urban surfaces.

# Fields
- `fRoofVeg`: Infiltration in first soil layer of vegetated roof [mm/h per horizontal roof area]
- `fGroundBare`: Infiltration in first soil layer of bare ground [mm/h per horizontal bare ground area]
- `fGroundVeg`: Infiltration in first soil layer of vegetated ground [mm/h per horizontal vegetated ground area]
- `fGroundImp`: Infiltration in impervious ground (usually zero) [mm/h per horizontal impervious ground area]
"""
Base.@kwdef mutable struct Infiltration{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    fRoofVeg::FT
    fGroundBare::FT
    fGroundVeg::FT
    fGroundImp::FT
end

function Infiltration(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Infiltration, Dict{String,Any}())
end

function roof_fields(::Type{Infiltration})
    return (:fRoofVeg,)
end

function ground_fields(::Type{Infiltration})
    return (:fGroundBare, :fGroundVeg, :fGroundImp)
end

"""
    Vwater{FT<:AbstractFloat, MR, MG} <: AbstractLayeredSoilVariables{FT}

Water volume in soil for different urban surfaces.

# Fields
- `VRoofSoilVeg`: Water volume in the different soil layers of roof [mm per horizontal roof area]
- `VGroundSoilImp`: Water volume in the different soil layers of ground under impervious [mm per horizontal impervious ground area]
- `VGroundSoilBare`: Water volume in the different soil layers of ground under bare [mm per horizontal bare ground area]
- `VGroundSoilVeg`: Water volume in the different soil layers of ground under vegetated [mm per horizontal vegetated ground area]
- `VGroundSoilTot`: Water volume in the different soil layers of ground total [mm per horizontal ground area]
"""
Base.@kwdef mutable struct Vwater{FT<:AbstractFloat,MR,MG} <:
                           AbstractLayeredSoilVariables{FT}
    VRoofSoilVeg::Vector{FT}
    VGroundSoilImp::Vector{FT}
    VGroundSoilBare::Vector{FT}
    VGroundSoilVeg::Vector{FT}
    VGroundSoilTot::Vector{FT}
end

function Vwater(::Type{FT}, soil::SoilParameters{FT}, args...) where {FT<:AbstractFloat}
    return initialize(
        FT, Vwater, Dict{String,Any}(), (FT, soil.roof.ms, soil.ground.ms), soil, args...
    )
end

function ground_fields(::Type{Vwater})
    return ["VGroundSoilImp", "VGroundSoilBare", "VGroundSoilVeg", "VGroundSoilTot"]
end

function roof_fields(::Type{Vwater})
    return ["VRoofSoilVeg"]
end

function Vwater(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    MR = length(data["VRoofSoilVeg"])
    MG = length(data["VGroundSoilImp"])
    return Vwater{FT,MR,MG}(;
        VRoofSoilVeg=data["VRoofSoilVeg"],
        VGroundSoilImp=data["VGroundSoilImp"],
        VGroundSoilBare=data["VGroundSoilBare"],
        VGroundSoilVeg=data["VGroundSoilVeg"],
        VGroundSoilTot=data["VGroundSoilTot"],
    )
end

function update!(
    dest::Vwater{FT,MR,MG}, src::Vwater{FT,MR,MG}
) where {FT<:AbstractFloat,MR,MG}
    dest.VRoofSoilVeg .= src.VRoofSoilVeg
    dest.VGroundSoilImp .= src.VGroundSoilImp
    dest.VGroundSoilBare .= src.VGroundSoilBare
    dest.VGroundSoilVeg .= src.VGroundSoilVeg
    dest.VGroundSoilTot .= src.VGroundSoilTot

    return nothing
end
"""
    dVwater_dt{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Change in water volume in soil for different urban surfaces.

# Fields
- `dVRoofSoilVeg_dt`: Change in water volume in the different soil layers of roof [mm per horizontal roof area]
- `dVGroundSoilImp_dt`: Change in water volume in the different soil layers of ground under impervious [mm per horizontal impervious ground area]
- `dVGroundSoilBare_dt`: Change in water volume in the different soil layers of ground under bare [mm per horizontal bare ground area]
- `dVGroundSoilVeg_dt`: Change in water volume in the different soil layers of ground under vegetated [mm per horizontal ground area]
- `dVGroundSoilTot_dt`: Change in water volume in the different soil layers of ground total [mm per horizontal ground area]
"""
# Same names as Vwater
Base.@kwdef mutable struct dVwater_dt{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    dVRoofSoilVeg_dt::FT
    dVGroundSoilImp_dt::FT
    dVGroundSoilBare_dt::FT
    dVGroundSoilVeg_dt::FT
    dVGroundSoilTot_dt::FT
end

function ground_fields(::Type{dVwater_dt})
    return (
        :dVGroundSoilImp_dt, :dVGroundSoilBare_dt, :dVGroundSoilVeg_dt, :dVGroundSoilTot_dt
    )
end

function roof_fields(::Type{dVwater_dt})
    return (:dVRoofSoilVeg_dt,)
end

function dVwater_dt(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, dVwater_dt, Dict{String,Any}())
end

"""
    Owater{FT<:AbstractFloat, MR, MG} <: AbstractLayeredSoilVariables{FT}

Soil moisture in different soil layers for urban surfaces.

# Fields
- `OwRoofSoilVeg`: Soil moisture in the different soil layers of roof [-]
- `OwGroundSoilImp`: Soil moisture in the different soil layers of ground under impervious [-]
- `OwGroundSoilBare`: Soil moisture in the different soil layers of ground under bare [-]
- `OwGroundSoilVeg`: Soil moisture in the different soil layers of ground under vegetated [-]
- `OwGroundSoilTot`: Soil moisture in the different soil layers of ground total [-]
"""
# Same names as Vwater
Base.@kwdef mutable struct Owater{FT<:AbstractFloat,MR,MG} <:
                           AbstractLayeredSoilVariables{FT}
    OwRoofSoilVeg::Vector{FT}
    OwGroundSoilImp::Vector{FT}
    OwGroundSoilBare::Vector{FT}
    OwGroundSoilVeg::Vector{FT}
    OwGroundSoilTot::Vector{FT}
end

function roof_fields(::Type{Owater})
    return ["OwRoofSoilVeg"]
end

function ground_fields(::Type{Owater})
    return ["OwGroundSoilImp", "OwGroundSoilBare", "OwGroundSoilVeg", "OwGroundSoilTot"]
end

function Owater(::Type{FT}, soil::SoilParameters{FT}, args...) where {FT<:AbstractFloat}
    return initialize(
        FT, Owater, Dict{String,Any}(), (FT, soil.roof.ms, soil.ground.ms), soil, args...
    )
end

function Owater(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    MR = length(data["OwRoofSoilVeg"])
    MG = length(data["OwGroundSoilImp"])
    return Owater{FT,MR,MG}(;
        OwRoofSoilVeg=data["OwRoofSoilVeg"],
        OwGroundSoilImp=data["OwGroundSoilImp"],
        OwGroundSoilBare=data["OwGroundSoilBare"],
        OwGroundSoilVeg=data["OwGroundSoilVeg"],
        OwGroundSoilTot=data["OwGroundSoilTot"],
    )
end

function update!(
    dest::Owater{FT,MR,MG}, src::Owater{FT,MR,MG}
) where {FT<:AbstractFloat,MR,MG}
    dest.OwRoofSoilVeg .= src.OwRoofSoilVeg
    dest.OwGroundSoilImp .= src.OwGroundSoilImp
    dest.OwGroundSoilBare .= src.OwGroundSoilBare
    dest.OwGroundSoilVeg .= src.OwGroundSoilVeg
    dest.OwGroundSoilTot .= src.OwGroundSoilTot

    return nothing
end

"""
    fix_soil_moisture!(
        dest::Owater{FT,MR,MG},
        roof::VegetatedSoilParameters{FT},
        ground::VegetatedSoilParameters{FT},
        O33::NamedTuple,
    ) where {FT<:AbstractFloat,MR,MG}

Fix soil moisture values in the soil layers according to the settings in the soil parameters.

# Arguments
- `dest::Owater{FT,MR,MG}`: Soil moisture variables to be modified
- `roof::VegetatedSoilParameters{FT}`: Roof soil parameters
- `ground::VegetatedSoilParameters{FT}`: Ground soil parameters
- `O33::NamedTuple`: Named tuple containing the fixed soil moisture values for roof and ground
"""
function fix_soil_moisture!(
    dest::Owater{FT,MR,MG},
    roof::VegetatedSoilParameters{FT},
    ground::VegetatedSoilParameters{FT},
    O33::NamedTuple,
) where {FT<:AbstractFloat,MR,MG}
    if roof.FixSM
        r = O33.roof
        SMReplace = fill(false, MR)
        SMReplace[roof.FixSM_LayerStart:roof.FixSM_LayerEnd] .= true
        dest.OwRoofSoilVeg[SMReplace .&& dest.OwRoofSoilVeg .< r] .= r
    end

    if ground.FixSM
        r = O33.ground
        SMReplace = fill(false, MG)
        SMReplace[ground.FixSM_LayerStart:ground.FixSM_LayerEnd] .= true
        dest.OwGroundSoilImp[SMReplace .&& dest.OwGroundSoilImp .< r] .= r
        dest.OwGroundSoilBare[SMReplace .&& dest.OwGroundSoilBare .< r] .= r
        dest.OwGroundSoilVeg[SMReplace .&& dest.OwGroundSoilVeg .< r] .= r
        dest.OwGroundSoilTot[SMReplace .&& dest.OwGroundSoilTot .< r] .= r
    end

    return nothing
end

"""
    OSwater{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Additional soil moisture variables for urban surfaces.

# Fields
- `OSwRoofSoilVeg`: Additional soil moisture values for roof soil layers [-]
- `OSwGroundSoilImp`: Additional soil moisture values for ground soil layers under impervious [-]
- `OSwGroundSoilBare`: Additional soil moisture values for ground soil layers under bare [-]
- `OSwGroundSoilVeg`: Additional soil moisture values for ground soil layers under vegetated [-]
- `OSwGroundSoilTot`: Additional soil moisture values for ground soil layers total [-]
"""
# Same names as Vwater
Base.@kwdef mutable struct OSwater{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    OSwRoofSoilVeg::FT
    OSwGroundSoilImp::FT
    OSwGroundSoilBare::FT
    OSwGroundSoilVeg::FT
    OSwGroundSoilTot::FT
end

function roof_fields(::Type{OSwater})
    return (:OSwRoofSoilVeg,)
end

function ground_fields(::Type{OSwater})
    return (:OSwGroundSoilImp, :OSwGroundSoilBare, :OSwGroundSoilVeg, :OSwGroundSoilTot)
end

function OSwater(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, OSwater, Dict{String,Any}())
end

"""
    Qinlat{FT<:AbstractFloat, MG} <: AbstractLayeredSoilVariables{FT}

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
Base.@kwdef mutable struct Qinlat{FT<:AbstractFloat,MG} <: AbstractLayeredSoilVariables{FT}
    Qin_bare2imp::Vector{FT}
    Qin_veg2imp::Vector{FT}
    Qin_veg2bare::Vector{FT}
    Qin_imp2bare::Vector{FT}
    Qin_bare2veg::Vector{FT}
    Qin_imp2veg::Vector{FT}
    Qin_imp::Vector{FT}
    Qin_bare::Vector{FT}
    Qin_veg::Vector{FT}
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

function Qinlat(::Type{FT}, soil::SoilParameters{FT}, args...) where {FT<:AbstractFloat}
    return initialize(FT, Qinlat, Dict{String,Any}(), (FT, soil.ground.ms), soil, args...)
end

function Qinlat(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    MG = length(data["Qin_bare2imp"])
    return Qinlat{FT,MG}(;
        Qin_bare2imp=data["Qin_bare2imp"],
        Qin_veg2imp=data["Qin_veg2imp"],
        Qin_veg2bare=data["Qin_veg2bare"],
        Qin_imp2bare=data["Qin_imp2bare"],
        Qin_bare2veg=data["Qin_bare2veg"],
        Qin_imp2veg=data["Qin_imp2veg"],
        Qin_imp=data["Qin_imp"],
        Qin_bare=data["Qin_bare"],
        Qin_veg=data["Qin_veg"],
    )
end

function update!(dest::Qinlat{FT,MG}, src::Qinlat{FT,MG}) where {FT<:AbstractFloat,MG}
    dest.Qin_bare2imp .= src.Qin_bare2imp
    dest.Qin_veg2imp .= src.Qin_veg2imp
    dest.Qin_veg2bare .= src.Qin_veg2bare
    dest.Qin_imp2bare .= src.Qin_imp2bare
    dest.Qin_bare2veg .= src.Qin_bare2veg
    dest.Qin_imp2veg .= src.Qin_imp2veg
    dest.Qin_imp .= src.Qin_imp
    dest.Qin_bare .= src.Qin_bare
    dest.Qin_veg .= src.Qin_veg

    return nothing
end

"""
    ExWater{FT<:AbstractFloat, MR, MG} <: AbstractLayeredSoilVariables{FT}

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
Base.@kwdef mutable struct ExWater{FT<:AbstractFloat,MR,MG} <:
                           AbstractLayeredSoilVariables{FT}
    ExWaterRoofVeg_H::Vector{FT}
    ExWaterRoofVeg_L::Vector{FT}
    ExWaterGroundImp_H::Vector{FT}
    ExWaterGroundImp_L::Vector{FT}
    ExWaterGroundBare_H::Vector{FT}
    ExWaterGroundBare_L::Vector{FT}
    ExWaterGroundVeg_H::Vector{FT}
    ExWaterGroundVeg_L::Vector{FT}
    ExWaterGroundTot_H::Vector{FT}
    ExWaterGroundTot_L::Vector{FT}
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

function ExWater(::Type{FT}, soil::SoilParameters{FT}, args...) where {FT<:AbstractFloat}
    return initialize(
        FT, ExWater, Dict{String,Any}(), (FT, soil.roof.ms, soil.ground.ms), soil, args...
    )
end

function ExWater(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    MR = length(data["ExWaterRoofVeg_H"])
    MG = length(data["ExWaterGroundImp_H"])
    return ExWater{FT,MR,MG}(;
        ExWaterRoofVeg_H=data["ExWaterRoofVeg_H"],
        ExWaterRoofVeg_L=data["ExWaterRoofVeg_L"],
        ExWaterGroundImp_H=data["ExWaterGroundImp_H"],
        ExWaterGroundImp_L=data["ExWaterGroundImp_L"],
        ExWaterGroundBare_H=data["ExWaterGroundBare_H"],
        ExWaterGroundBare_L=data["ExWaterGroundBare_L"],
        ExWaterGroundVeg_H=data["ExWaterGroundVeg_H"],
        ExWaterGroundVeg_L=data["ExWaterGroundVeg_L"],
        ExWaterGroundTot_H=data["ExWaterGroundTot_H"],
        ExWaterGroundTot_L=data["ExWaterGroundTot_L"],
    )
end

function update!(
    dest::ExWater{FT,MR,MG}, src::ExWater{FT,MR,MG}
) where {FT<:AbstractFloat,MR,MG}
    dest.ExWaterRoofVeg_H .= src.ExWaterRoofVeg_H
    dest.ExWaterRoofVeg_L .= src.ExWaterRoofVeg_L
    dest.ExWaterGroundImp_H .= src.ExWaterGroundImp_H
    dest.ExWaterGroundImp_L .= src.ExWaterGroundImp_L
    dest.ExWaterGroundBare_H .= src.ExWaterGroundBare_H
    dest.ExWaterGroundBare_L .= src.ExWaterGroundBare_L
    dest.ExWaterGroundVeg_H .= src.ExWaterGroundVeg_H
    dest.ExWaterGroundVeg_L .= src.ExWaterGroundVeg_L
    dest.ExWaterGroundTot_H .= src.ExWaterGroundTot_H
    dest.ExWaterGroundTot_L .= src.ExWaterGroundTot_L

    return nothing
end

"""
    SoilPotW{FT<:AbstractFloat} <: AbstractModelVariables{FT}

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
Base.@kwdef mutable struct SoilPotW{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    SoilPotWRoofVeg_H::FT
    SoilPotWRoofVeg_L::FT
    SoilPotWGroundImp_H::FT
    SoilPotWGroundImp_L::FT
    SoilPotWGroundBare_H::FT
    SoilPotWGroundBare_L::FT
    SoilPotWGroundVeg_H::FT
    SoilPotWGroundVeg_L::FT
    SoilPotWGroundTot_H::FT
    SoilPotWGroundTot_L::FT
end

function SoilPotW(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, SoilPotW, Dict{String,Any}())
end

function SoilPotW(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return SoilPotW{FT}(;
        SoilPotWRoofVeg_H=data["SoilPotWRoofVeg_H"],
        SoilPotWRoofVeg_L=data["SoilPotWRoofVeg_L"],
        SoilPotWGroundImp_H=data["SoilPotWGroundImp_H"],
        SoilPotWGroundImp_L=data["SoilPotWGroundImp_L"],
        SoilPotWGroundBare_H=data["SoilPotWGroundBare_H"],
        SoilPotWGroundBare_L=data["SoilPotWGroundBare_L"],
        SoilPotWGroundVeg_H=data["SoilPotWGroundVeg_H"],
        SoilPotWGroundVeg_L=data["SoilPotWGroundVeg_L"],
        SoilPotWGroundTot_H=data["SoilPotWGroundTot_H"],
        SoilPotWGroundTot_L=data["SoilPotWGroundTot_L"],
    )
end

function update!(dest::SoilPotW{FT}, src::SoilPotW{FT}) where {FT<:AbstractFloat}
    dest.SoilPotWRoofVeg_H = src.SoilPotWRoofVeg_H
    dest.SoilPotWRoofVeg_L = src.SoilPotWRoofVeg_L
    dest.SoilPotWGroundImp_H = src.SoilPotWGroundImp_H
    dest.SoilPotWGroundImp_L = src.SoilPotWGroundImp_L
    dest.SoilPotWGroundBare_H = src.SoilPotWGroundBare_H
    dest.SoilPotWGroundBare_L = src.SoilPotWGroundBare_L
    dest.SoilPotWGroundVeg_H = src.SoilPotWGroundVeg_H
    dest.SoilPotWGroundVeg_L = src.SoilPotWGroundVeg_L
    dest.SoilPotWGroundTot_H = src.SoilPotWGroundTot_H
    dest.SoilPotWGroundTot_L = src.SoilPotWGroundTot_L

    return nothing
end

function ground_fields(::Type{SoilPotW})
    return (
        :SoilPotWGroundImp_H,
        :SoilPotWGroundImp_L,
        :SoilPotWGroundBare_H,
        :SoilPotWGroundBare_L,
        :SoilPotWGroundVeg_H,
        :SoilPotWGroundVeg_L,
        :SoilPotWGroundTot_H,
        :SoilPotWGroundTot_L,
    )
end

function roof_fields(::Type{SoilPotW})
    return (:SoilPotWRoofVeg_H, :SoilPotWRoofVeg_L)
end

"""
    CiCO2Leaf{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Intercellular CO2 concentration in leaf for different urban surfaces.

# Fields
- `CiCO2LeafRoofVegSun`: Ci_sun_veg sunlit roof leafs [umolCO2/mol]
- `CiCO2LeafRoofVegShd`: Ci_shd_veg shaded roof leafs [umolCO2/mol]
- `CiCO2LeafGroundVegSun`: Ci_sun_veg sunlit ground leafs [umolCO2/mol]
- `CiCO2LeafGroundVegShd`: Ci_shd_veg shaded ground leafs [umolCO2/mol]
- `CiCO2LeafTreeSun`: Ci_sun_veg sunlit tree leafs [umolCO2/mol]
- `CiCO2LeafTreeShd`: Ci_shd_veg shaded tree leafs [umolCO2/mol]
"""
Base.@kwdef mutable struct CiCO2Leaf{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    CiCO2LeafRoofVegSun::FT
    CiCO2LeafRoofVegShd::FT
    CiCO2LeafGroundVegSun::FT
    CiCO2LeafGroundVegShd::FT
    CiCO2LeafTreeSun::FT
    CiCO2LeafTreeShd::FT
end

function CiCO2Leaf(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, CiCO2Leaf, Dict{String,Any}(), (FT,))
end

function CiCO2Leaf(::Type{FT}, data::AbstractDict) where {FT<:AbstractFloat}
    return CiCO2Leaf{FT}(;
        CiCO2LeafRoofVegSun=data["CiCO2LeafRoofVegSun"],
        CiCO2LeafRoofVegShd=data["CiCO2LeafRoofVegShd"],
        CiCO2LeafGroundVegSun=data["CiCO2LeafGroundVegSun"],
        CiCO2LeafGroundVegShd=data["CiCO2LeafGroundVegShd"],
        CiCO2LeafTreeSun=data["CiCO2LeafTreeSun"],
        CiCO2LeafTreeShd=data["CiCO2LeafTreeShd"],
    )
end

function update!(dest::CiCO2Leaf{FT}, src::CiCO2Leaf{FT}) where {FT<:AbstractFloat}
    dest.CiCO2LeafRoofVegSun = src.CiCO2LeafRoofVegSun
    dest.CiCO2LeafRoofVegShd = src.CiCO2LeafRoofVegShd
    dest.CiCO2LeafGroundVegSun = src.CiCO2LeafGroundVegSun
    dest.CiCO2LeafGroundVegShd = src.CiCO2LeafGroundVegShd
    dest.CiCO2LeafTreeSun = src.CiCO2LeafTreeSun
    dest.CiCO2LeafTreeShd = src.CiCO2LeafTreeShd

    return nothing
end

function roof_fields(::Type{CiCO2Leaf})
    return (:CiCO2LeafRoofVegSun, :CiCO2LeafRoofVegShd)
end

function ground_fields(::Type{CiCO2Leaf})
    return (
        :CiCO2LeafGroundVegSun, :CiCO2LeafGroundVegShd, :CiCO2LeafTreeSun, :CiCO2LeafTreeShd
    )
end
"""
    WaterFluxVariables{FT<:AbstractFloat, MR, MG} <: AbstractModelVariableSet{FT}

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
Base.@kwdef struct WaterFluxVariables{FT<:AbstractFloat,MR,MG} <:
                   AbstractModelVariableSet{FT}
    Eflux::Eflux{FT}
    Runoff::Runoff{FT}
    Runon::Runon{FT}
    Leakage::Leakage{FT}
    Interception::Interception{FT}
    dInt_dt::dInt_dt{FT}
    Infiltration::Infiltration{FT}
    Vwater::Vwater{FT,MR,MG}
    dVwater_dt::dVwater_dt{FT}
    Owater::Owater{FT,MR,MG}
    OSwater::OSwater{FT}
    Qinlat::Qinlat{FT,MG}
    ExWater::ExWater{FT,MR,MG}
    SoilPotW::SoilPotW{FT}
    CiCO2Leaf::CiCO2Leaf{FT}
end

function WaterFluxVariables(
    ::Type{FT}, soil_parameters::SoilParameters{FT}
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        WaterFluxVariables,
        Dict{String,Any}(),
        (FT, soil_parameters.roof.ms, soil_parameters.ground.ms),
        soil_parameters,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{WaterFluxVariables},
    data::Dict{String,Any},
    params::Tuple,
    soil::SoilParameters{FT},
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["Eflux"] = Eflux(FT)
    processed["Runoff"] = Runoff(FT)
    processed["Runon"] = Runon(FT)
    processed["Leakage"] = Leakage(FT)
    processed["Interception"] = Interception(FT)
    processed["dInt_dt"] = dInt_dt(FT)
    processed["Infiltration"] = Infiltration(FT)
    processed["Vwater"] = Vwater(FT, soil)
    processed["dVwater_dt"] = dVwater_dt(FT)
    processed["Owater"] = Owater(FT, soil)
    processed["OSwater"] = OSwater(FT)
    processed["Qinlat"] = Qinlat(FT, soil)
    processed["ExWater"] = ExWater(FT, soil)
    processed["SoilPotW"] = SoilPotW(FT)
    processed["CiCO2Leaf"] = CiCO2Leaf(FT)

    return processed
end

function ModelComponents.outputs_to_save(::Type{WaterFluxVariables}, ::Type{PlotOutputs})
    return (:Runoff, :Runon, :Leakage, :Interception, :dInt_dt, :dVwater_dt, :Owater)
end

function ModelComponents.outputs_to_save(
    ::Type{WaterFluxVariables}, ::Type{ExtendedEnergyClimateOutputs}
)
    return (:Eflux, :Infiltration, :Vwater, :SoilPotW)
end

function ModelComponents.outputs_to_save(
    ::Type{WaterFluxVariables}, ::Type{ExtendedOutputs}
)
    return (:OSwater, :Qinlat, :ExWater, :CiCO2Leaf)
end

function update!(
    x::WaterFluxVariables{FT}, results::NamedTuple, fn::EBWBRoofDispatcher
) where {FT<:AbstractFloat}
    for field in fieldnames(WaterFluxVariables{FT})
        field_value = getfield(x, field)
        field_type = typeof(field_value).name.wrapper
        _update!(field_value, results, roof_fields(field_type))
    end

    return nothing
end

function update!(
    x::WaterFluxVariables{FT}, results::NamedTuple, ::EBWBCanyonDispatcher
) where {FT<:AbstractFloat}
    for field in fieldnames(WaterFluxVariables{FT})
        field_value = getfield(x, field)
        field_type = typeof(field_value).name.wrapper
        _update!(field_value, results, ground_fields(field_type))
    end

    return nothing
end
