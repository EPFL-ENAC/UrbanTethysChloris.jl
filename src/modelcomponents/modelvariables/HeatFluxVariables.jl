"""
    Hflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Sensible heat fluxes for different urban surfaces.

# Fields
- `HfluxRoofImp`: Sensible heat flux of impervious roof area to atmosphere [W/m² horizontal impervious roof area]
- `HfluxRoofVeg`: Sensible heat flux of vegetated roof area to atmosphere [W/m² horizontal vegetated roof area]
- `HfluxRoof`: Total sensible heat flux of roof area to atmosphere [W/m² horizontal roof area]
- `HfluxGroundImp`: Sensible heat flux of impervious ground area to canyon [W/m² horizontal impervious ground area]
- `HfluxGroundBare`: Sensible heat flux of bare ground area to canyon [W/m² horizontal bare ground area]
- `HfluxGroundVeg`: Sensible heat flux of vegetated ground area to canyon [W/m² horizontal vegetated ground area]
- `HfluxGround`: Sensible heat flux of ground area to canyon [W/m² horizontal ground area]
- `HfluxTree`: Sensible heat flux of tree canopy to canyon [W/m² horizontally projected tree area]
- `HfluxWallSun`: Sensible heat flux of sunlit wall to canyon [W/m² vertical wall area]
- `HfluxWallShade`: Sensible heat flux of shaded wall to canyon [W/m² vertical wall area]
- `HfluxCanyon`: Sensible heat flux of canyon to atmosphere [W/m² horizontal canyon area]
- `HfluxUrban`: Total sensible heat flux of urban area to atmosphere [W/m² horizontal urban area]
- `dS_H_air`: Change in sensible heat storage in canyon air volume due to temperature change [W/m² horizontal area]
"""
Base.@kwdef mutable struct Hflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    HfluxRoofImp::FT
    HfluxRoofVeg::FT
    HfluxRoof::FT
    HfluxGroundImp::FT
    HfluxGroundBare::FT
    HfluxGroundVeg::FT
    HfluxGround::FT
    HfluxTree::FT
    HfluxWallSun::FT
    HfluxWallShade::FT
    HfluxCanyon::FT
    HfluxUrban::FT
    dS_H_air::FT
end

function Hflux(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Hflux, Dict{String,Any}())
end

"""
    LEflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Latent heat fluxes for different urban surfaces.

# Fields
- `LEfluxRoofImp`: Latent heat flux of intercepted water from impervious roof area to atmosphere [W/m² horizontal impervious roof area]
- `LEfluxRoofVegInt`: Latent heat flux of intercepted water on roof vegetation to atmosphere [W/m² horizontal vegetated roof area]
- `LEfluxRoofVegPond`: Latent heat flux of intercepted water on ground under roof vegetation to atmosphere [W/m² horizontal vegetated roof area]
- `LEfluxRoofVegSoil`: Latent heat flux of water from roof soil under vegetation to atmosphere [W/m² horizontal vegetated roof area]
- `LTEfluxRoofVeg`: Latent heat flux of transpiration from roof plants to atmosphere [W/m² horizontal vegetated roof area]
- `LEfluxRoofVeg`: Total latent heat flux of vegetated roof to atmosphere [W/m² horizontal vegetated roof area]
- `LEfluxRoof`: Total latent heat flux of roof to atmosphere [W/m² horizontal roof area]
- `LEfluxGroundImp`: Latent heat flux of intercepted water on impervious ground area to canyon [W/m² horizontal impervious ground area]
- `LEfluxGroundBarePond`: Latent heat flux of water on bare ground to canyon [W/m² horizontal bare ground area]
- `LEfluxGroundBareSoil`: Latent heat flux of water from bare ground to canyon [W/m² horizontal bare ground area]
- `LEfluxGroundBare`: Total latent heat flux of bare ground area to canyon [W/m² horizontal bare ground area]
- `LEfluxGroundVegInt`: Latent heat flux of intercepted water on ground vegetation to canyon [W/m² horizontal vegetated ground area]
- `LEfluxGroundVegPond`: Latent heat flux of intercepted water on ground under vegetation to canyon [W/m² horizontal vegetated ground area]
- `LEfluxGroundVegSoil`: Latent heat flux of water from ground soil under vegetation to canyon [W/m² horizontal vegetated ground area]
- `LTEfluxGroundVeg`: Latent heat flux of transpiration from ground plants to canyon [W/m² horizontal vegetated ground area]
- `LEfluxGroundVeg`: Total latent heat flux of vegetated ground to canyon [W/m² horizontal vegetated ground area]
- `LEfluxGround`: Total latent heat flux of ground to canyon [W/m² horizontal ground area]
- `LEfluxTreeInt`: Latent heat flux of intercepted water on tree canopy to canyon [W/m² horizontally projected tree area]
- `LTEfluxTree`: Latent heat flux of transpiration from tree canopy to canyon [W/m² horizontally projected tree area]
- `LEfluxTree`: Total latent heat flux of tree canopy to canyon [W/m² horizontally projected tree area]
- `LEfluxWallSun`: Latent heat flux of sunlit wall to canyon [W/m² vertical wall area]
- `LEfluxWallShade`: Latent heat flux of shaded wall to canyon [W/m² vertical wall area]
- `LEfluxCanyon`: Latent heat flux of canyon to atmosphere [W/m² horizontal canyon area]
- `LEfluxUrban`: Total latent heat flux of urban area to atmosphere [W/m² horizontal urban area]
- `dS_LE_air`: Change in latent heat storage in canyon air volume due to moisture change [W/m² horizontal area]
"""
Base.@kwdef mutable struct LEflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    LEfluxRoofImp::FT
    LEfluxRoofVegInt::FT
    LEfluxRoofVegPond::FT
    LEfluxRoofVegSoil::FT
    LTEfluxRoofVeg::FT
    LEfluxRoofVeg::FT
    LEfluxRoof::FT
    LEfluxGroundImp::FT
    LEfluxGroundBarePond::FT
    LEfluxGroundBareSoil::FT
    LEfluxGroundBare::FT
    LEfluxGroundVegInt::FT
    LEfluxGroundVegPond::FT
    LEfluxGroundVegSoil::FT
    LTEfluxGroundVeg::FT
    LEfluxGroundVeg::FT
    LEfluxGround::FT
    LEfluxTreeInt::FT
    LTEfluxTree::FT
    LEfluxTree::FT
    LEfluxWallSun::FT
    LEfluxWallShade::FT
    LEfluxCanyon::FT
    LEfluxUrban::FT
    dS_LE_air::FT
end

function LEflux(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, LEflux, Dict{String,Any}())
end

"""
    Gflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Conductive heat fluxes for different urban surfaces.

# Fields
- `G1RoofImp`: Conductive heat flux of first layer of impervious roof [W/m² horizontal impervious roof area]
- `G1RoofVeg`: Conductive heat flux of first layer of vegetated roof [W/m² horizontal vegetated roof area]
- `G2RoofImp`: Conductive heat flux of second layer of impervious roof [W/m² horizontal impervious roof area]
- `G2RoofVeg`: Conductive heat flux of second layer of vegetated roof [W/m² horizontal vegetated roof area]
- `G1Roof`: Total conductive heat flux of first layer of roof [W/m² horizontal roof area]
- `G2Roof`: Total conductive heat flux of second layer of roof [W/m² horizontal roof area]
- `G1GroundImp`: Conductive heat flux of impervious ground [W/m² horizontal impervious ground area]
- `G1GroundBare`: Conductive heat flux of bare ground [W/m² horizontal bare ground area]
- `G1GroundVeg`: Conductive heat flux of vegetated ground [W/m² horizontal vegetated ground area]
- `G1Ground`: Total conductive heat flux of ground [W/m² horizontal ground area]
- `GTree`: Conductive heat flux tree [W/m² horizontally projected tree area]
- `G1WallSun`: Conductive heat flux of first layer of sunlit wall [W/m² vertical wall area]
- `G1WallShade`: Conductive heat flux of first layer of shaded wall [W/m² vertical wall area]
- `G2WallSun`: Conductive heat flux of second layer of sunlit wall [W/m² vertical wall area]
- `G2WallShade`: Conductive heat flux of second layer of shaded wall [W/m² vertical wall area]
- `G1Canyon`: Total conductive heat flux G1 (walls and ground) area-averaged per m² canyon ground [W/m² horizontal canyon area]
- `G2Canyon`: Total conductive heat flux G2 (walls and ground=0) area-averaged per m² canyon ground [W/m² horizontal canyon area]
- `G1Urban`: Total conductive heat flux G1 area-averaged per m² urban [W/m² horizontal urban area]
- `G2Urban`: Total conductive heat flux G2 (G2 of ground is 0) area-averaged per m² urban [W/m² horizontal urban area]
"""
Base.@kwdef mutable struct Gflux{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    G1RoofImp::FT
    G1RoofVeg::FT
    G2RoofImp::FT
    G2RoofVeg::FT
    G1Roof::FT
    G2Roof::FT
    G1GroundImp::FT
    G1GroundBare::FT
    G1GroundVeg::FT
    G1Ground::FT
    GTree::FT
    G1WallSun::FT
    G1WallShade::FT
    G2WallSun::FT
    G2WallShade::FT
    G1Canyon::FT
    G2Canyon::FT
    G1Urban::FT
    G2Urban::FT
end

function Gflux(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Gflux, Dict{String,Any}())
end

"""
    dStorage{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Heat storage in different urban surfaces.

# Fields
- `dsRoofImp`: Storage of energy in impervious roof [W/m² horizontal impervious roof area]
- `dsRoofVeg`: Storage of energy in vegetated roof [W/m² horizontal vegetated roof area]
- `dsRoof`: Storage of energy in total roof [W/m² horizontal roof area]
- `dsGroundImp`: Storage of energy in impervious ground [W/m² horizontal impervious ground area]
- `dsGroundBare`: Storage of energy in bare ground [W/m² horizontal bare ground area]
- `dsGroundVeg`: Storage of energy in vegetated ground [W/m² horizontal vegetated ground area]
- `dsTree`: Storage of energy in tree canopy [W/m² horizontally projected tree area]
- `dsWallSun`: Storage of energy in sunlit wall [W/m² vertical wall area]
- `dsWallShade`: Storage of energy in shaded wall [W/m² vertical wall area]
- `dsCanyonAir`: Storage of energy in canyon air [W/m² horizontal canyon area]
"""
Base.@kwdef mutable struct dStorage{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    dsRoofImp::FT
    dsRoofVeg::FT
    dsRoof::FT
    dsGroundImp::FT
    dsGroundBare::FT
    dsGroundVeg::FT
    dsTree::FT
    dsWallSun::FT
    dsWallShade::FT
    dsCanyonAir::FT
end

function dStorage(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, dStorage, Dict{String,Any}())
end

"""
    Results2mEnergyFlux{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Energy fluxes at 2m canyon height.

# Fields
- `DHi`: Sensible heat flux from impervious surfaces at 2m height [W/m²]
- `Himp_2m`: Sensible heat flux from impervious ground at 2m height [W/m²]
- `Hbare_2m`: Sensible heat flux from bare ground at 2m height [W/m²]
- `Hveg_2m`: Sensible heat flux from vegetated ground at 2m height [W/m²]
- `Hwsun_2m`: Sensible heat flux from sunlit wall at 2m height [W/m²]
- `Hwshade_2m`: Sensible heat flux from shaded wall at 2m height [W/m²]
- `Hcan_2m`: Total sensible heat flux at 2m height [W/m²]
- `DEi`: Latent heat flux at 2m height [W/m²]
- `Eimp_2m`: Latent heat flux from impervious ground at 2m height [W/m²]
- `Ebare_soil_2m`: Latent heat flux from bare soil at 2m height [W/m²]
- `Eveg_int_2m`: Latent heat flux from vegetation interception at 2m height [W/m²]
- `Eveg_soil_2m`: Latent heat flux from vegetated soil at 2m height [W/m²]
- `TEveg_2m`: Transpiration flux at 2m height [W/m²]
- `Ecan_2m`: Total latent heat flux at 2m height [W/m²]
"""
Base.@kwdef mutable struct Results2mEnergyFlux{FT<:AbstractFloat} <:
                           AbstractModelVariables{FT}
    DHi::FT
    Himp_2m::FT
    Hbare_2m::FT
    Hveg_2m::FT
    Hwsun_2m::FT
    Hwshade_2m::FT
    Hcan_2m::FT
    DEi::FT
    Eimp_2m::FT
    Ebare_soil_2m::FT
    Eveg_int_2m::FT
    Eveg_soil_2m::FT
    TEveg_2m::FT
    Ecan_2m::FT
end

function Results2mEnergyFlux(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, Results2mEnergyFlux, Dict{String,Any}())
end

"""
    HeatFluxVariables{FT<:AbstractFloat} <: AbstractModelVariablesSet{FT}

Container for all heat flux variable components.

# Fields
- `Hflux`: Sensible heat fluxes for different urban surfaces
- `LEflux`: Latent heat fluxes for different urban surfaces
- `Gflux`: Conductive heat fluxes for different urban surfaces
- `dStorage`: Heat storage in different urban surfaces
- `Results2mEnergyFlux`: Energy fluxes at 2m canyon height
"""
Base.@kwdef struct HeatFluxVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}
    Hflux::Hflux{FT}
    LEflux::LEflux{FT}
    Gflux::Gflux{FT}
    dStorage::dStorage{FT}
    Results2mEnergyFlux::Results2mEnergyFlux{FT}
end

function HeatFluxVariables(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, HeatFluxVariables, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{HeatFluxVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()

    processed["Hflux"] = Hflux(FT)
    processed["LEflux"] = LEflux(FT)
    processed["Gflux"] = Gflux(FT)
    processed["dStorage"] = dStorage(FT)
    processed["Results2mEnergyFlux"] = Results2mEnergyFlux(FT)

    return processed
end

function ModelComponents.outputs_to_save(
    ::Type{HeatFluxVariables}, ::Type{ExtendedEnergyClimateOutputs}
)
    return (:Hflux, :LEflux, :Gflux, :dStorage, :Results2mEnergyFlux)
end
