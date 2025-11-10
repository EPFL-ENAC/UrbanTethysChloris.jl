"""
    Hflux{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct Hflux{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    HfluxRoofImp::Array{FT,N}
    HfluxRoofVeg::Array{FT,N}
    HfluxRoof::Array{FT,N}
    HfluxGroundImp::Array{FT,N}
    HfluxGroundBare::Array{FT,N}
    HfluxGroundVeg::Array{FT,N}
    HfluxGround::Array{FT,N}
    HfluxTree::Array{FT,N}
    HfluxWallSun::Array{FT,N}
    HfluxWallShade::Array{FT,N}
    HfluxCanyon::Array{FT,N}
    HfluxUrban::Array{FT,N}
    dS_H_air::Array{FT,N}
end

function initialize_hflux(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Hflux, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_hflux(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, Hflux, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    LEflux{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct LEflux{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    LEfluxRoofImp::Array{FT,N}
    LEfluxRoofVegInt::Array{FT,N}
    LEfluxRoofVegPond::Array{FT,N}
    LEfluxRoofVegSoil::Array{FT,N}
    LTEfluxRoofVeg::Array{FT,N}
    LEfluxRoofVeg::Array{FT,N}
    LEfluxRoof::Array{FT,N}
    LEfluxGroundImp::Array{FT,N}
    LEfluxGroundBarePond::Array{FT,N}
    LEfluxGroundBareSoil::Array{FT,N}
    LEfluxGroundBare::Array{FT,N}
    LEfluxGroundVegInt::Array{FT,N}
    LEfluxGroundVegPond::Array{FT,N}
    LEfluxGroundVegSoil::Array{FT,N}
    LTEfluxGroundVeg::Array{FT,N}
    LEfluxGroundVeg::Array{FT,N}
    LEfluxGround::Array{FT,N}
    LEfluxTreeInt::Array{FT,N}
    LTEfluxTree::Array{FT,N}
    LEfluxTree::Array{FT,N}
    LEfluxWallSun::Array{FT,N}
    LEfluxWallShade::Array{FT,N}
    LEfluxCanyon::Array{FT,N}
    LEfluxUrban::Array{FT,N}
    dS_LE_air::Array{FT,N}
end

function initialize_leflux(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, LEflux, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_leflux(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, LEflux, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    Gflux{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct Gflux{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    G1RoofImp::Array{FT,N}
    G1RoofVeg::Array{FT,N}
    G2RoofImp::Array{FT,N}
    G2RoofVeg::Array{FT,N}
    G1Roof::Array{FT,N}
    G2Roof::Array{FT,N}
    G1GroundImp::Array{FT,N}
    G1GroundBare::Array{FT,N}
    G1GroundVeg::Array{FT,N}
    G1Ground::Array{FT,N}
    GTree::Array{FT,N}
    G1WallSun::Array{FT,N}
    G1WallShade::Array{FT,N}
    G2WallSun::Array{FT,N}
    G2WallShade::Array{FT,N}
    G1Canyon::Array{FT,N}
    G2Canyon::Array{FT,N}
    G1Urban::Array{FT,N}
    G2Urban::Array{FT,N}
end

function initialize_gflux(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, Gflux, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_gflux(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, Gflux, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    dStorage{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct dStorage{FT<:AbstractFloat,N} <: Abstract1PModelVariables{FT,N}
    dsRoofImp::Array{FT,N}
    dsRoofVeg::Array{FT,N}
    dsRoof::Array{FT,N}
    dsGroundImp::Array{FT,N}
    dsGroundBare::Array{FT,N}
    dsGroundVeg::Array{FT,N}
    dsTree::Array{FT,N}
    dsWallSun::Array{FT,N}
    dsWallShade::Array{FT,N}
    dsCanyonAir::Array{FT,N}
end

function initialize_dstorage(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(FT, dStorage, Dict{String,Any}(), (FT, dimension_value(TimeSlice())))
end

function initialize_dstorage(::Type{FT}, ::TimeSeries, hours::Int) where {FT<:AbstractFloat}
    return initialize(
        FT, dStorage, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    Results2mEnergyFlux{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct Results2mEnergyFlux{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    DHi::Array{FT,N}
    Himp_2m::Array{FT,N}
    Hbare_2m::Array{FT,N}
    Hveg_2m::Array{FT,N}
    Hwsun_2m::Array{FT,N}
    Hwshade_2m::Array{FT,N}
    Hcan_2m::Array{FT,N}
    DEi::Array{FT,N}
    Eimp_2m::Array{FT,N}
    Ebare_soil_2m::Array{FT,N}
    Eveg_int_2m::Array{FT,N}
    Eveg_soil_2m::Array{FT,N}
    TEveg_2m::Array{FT,N}
    Ecan_2m::Array{FT,N}
end

function initialize_results2m_energy_flux(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, Results2mEnergyFlux, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_results2m_energy_flux(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        Results2mEnergyFlux,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
    )
end

"""
    HeatFluxVariables{FT<:AbstractFloat, N} <: Abstract1PModelVariablesSet{FT, N}

Container for all heat flux variable components.

# Fields
- `Hflux`: Sensible heat fluxes for different urban surfaces
- `LEflux`: Latent heat fluxes for different urban surfaces
- `Gflux`: Conductive heat fluxes for different urban surfaces
- `dStorage`: Heat storage in different urban surfaces
- `Results2mEnergyFlux`: Energy fluxes at 2m canyon height
"""
Base.@kwdef struct HeatFluxVariables{FT<:AbstractFloat,N} <:
                   Abstract1PModelVariablesSet{FT,N}
    Hflux::Hflux{FT,N}
    LEflux::LEflux{FT,N}
    Gflux::Gflux{FT,N}
    dStorage::dStorage{FT,N}
    Results2mEnergyFlux::Results2mEnergyFlux{FT,N}
end

function initialize_heat_flux_variables(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, HeatFluxVariables, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{HeatFluxVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["Hflux"] = initialize_hflux(FT, dimensionality_type(params[2]))
    processed["LEflux"] = initialize_leflux(FT, dimensionality_type(params[2]))
    processed["Gflux"] = initialize_gflux(FT, dimensionality_type(params[2]))
    processed["dStorage"] = initialize_dstorage(FT, dimensionality_type(params[2]))
    processed["Results2mEnergyFlux"] = initialize_results2m_energy_flux(
        FT, dimensionality_type(params[2])
    )
    return processed
end

function initialize_heat_flux_variables(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        HeatFluxVariables,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{HeatFluxVariables}, data::Dict{String,Any}, params::Tuple, hours::Int
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["Hflux"] = initialize_hflux(FT, dimensionality_type(params[2]), hours)
    processed["LEflux"] = initialize_leflux(FT, dimensionality_type(params[2]), hours)
    processed["Gflux"] = initialize_gflux(FT, dimensionality_type(params[2]), hours)
    processed["dStorage"] = initialize_dstorage(FT, dimensionality_type(params[2]), hours)
    processed["Results2mEnergyFlux"] = initialize_results2m_energy_flux(
        FT, dimensionality_type(params[2]), hours
    )
    return processed
end
