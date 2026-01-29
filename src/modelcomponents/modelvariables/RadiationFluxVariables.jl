"""
    AbsorbedRadiationFluxVariablesSubset{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Absorbed shortwave radiation for different urban surfaces.

# Fields
- `RoofImp`: Absorbed shortwave radiation roof impervious area [W/m² horizontal impervious roof area]
- `RoofVeg`: Absorbed shortwave radiation roof vegetated area [W/m² horizontal vegetated roof area]
- `TotalRoof`: Total absorbed shortwave radiation by the roof area [W/m² horizontal roof area]
- `GroundImp`: Absorbed shortwave radiation ground impervious area [W/m² horizontal impervious ground area]
- `GroundBare`: Absorbed shortwave radiation ground bare area [W/m² horizontal bare ground area]
- `GroundVeg`: Absorbed shortwave radiation ground vegetated area [W/m² horizontal vegetated ground area]
- `Tree`: Absorbed shortwave radiation tree canopy [W/m² horizontally projected tree area]
- `WallSun`: Absorbed shortwave radiation sunlit area [W/m² vertical wall area]
- `WallShade`: Absorbed shortwave radiation shaded area [W/m² vertical wall area]
- `TotalGround`: Total absorbed shortwave radiation by the canyon ground area [W/m²]
- `TotalCanyon`: Total absorbed shortwave radiation by all the canyon facets [W/m²]
- `TotalUrban`: Total absorbed shortwave radiation by all the urban elements (roof plus canyon) [W/m²]
- `WallSunExt`: Absorbed shortwave radiation by the exterior sunlit wall area, average over wall and window area [W/m² vertical wall area]
- `WallShadeExt`: Absorbed shortwave radiation by the exterior sunlit wall area, average over wall and window area [W/m² vertical wall area]
- `WallSunTransmitted`: Shortwave radiation transmitted through the windows on the sunlit wall area [W/m² vertical wall area]
- `WallShadeTransmitted`: Shortwave radiation transmitted through the windows on the sunlit wall area [W/m² vertical wall area]
"""
Base.@kwdef mutable struct AbsorbedRadiationFluxVariablesSubset{FT<:AbstractFloat} <:
                           AbstractModelVariables{FT}
    SWRabsRoofImp::FT
    RoofVeg::FT
    TotalRoof::FT
    GroundImp::FT
    GroundBare::FT
    GroundVeg::FT
    Tree::FT
    WallSun::FT
    WallShade::FT
    TotalGround::FT
    TotalCanyon::FT
    TotalUrban::FT
    WallSunExt::FT
    WallShadeExt::FT
    WallSunTransmitted::FT
    WallShadeTransmitted::FT
end

function AbsorbedRadiationFluxVariablesSubset(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, AbsorbedRadiationFluxVariablesSubset, Dict{String,Any}())
end

"""
    DefaultRadiationFluxVariablesSubset{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Radiation flux for different urban surfaces.

# Fields
- `RoofImp`: Incoming shortwave radiation roof impervious area [W/m² horizontal impervious roof area]
- `RoofVeg`: Incoming shortwave radiation roof vegetated area [W/m² horizontal vegetated roof area]
- `TotalRoof`: Total incoming shortwave radiation by the roof area [W/m² horizontal roof area]
- `GroundImp`: Incoming shortwave radiation ground impervious area [W/m² horizontal impervious ground area]
- `GroundBare`: Incoming shortwave radiation ground bare area [W/m² horizontal bare ground area]
- `GroundVeg`: Incoming shortwave radiation ground vegetated area [W/m² horizontal vegetated ground area]
- `Tree`: Incoming shortwave radiation tree canopy [W/m² horizontally projected tree area]
- `WallSun`: Incoming shortwave radiation sunlit area [W/m² vertical wall area]
- `WallShade`: Incoming shortwave radiation shaded area [W/m² vertical wall area]
- `TotalGround`: Total incoming shortwave radiation by the canyon ground area [W/m²]
- `TotalCanyon`: Total incoming shortwave radiation by all the canyon facets [W/m²]
- `TotalUrban`: Total incoming shortwave radiation by all the urban elements (roof plus canyon) [W/m²]
"""
Base.@kwdef mutable struct DefaultRadiationFluxVariablesSubset{FT<:AbstractFloat} <:
                           AbstractModelVariables{FT}
    RoofImp::FT
    RoofVeg::FT
    TotalRoof::FT
    GroundImp::FT
    GroundBare::FT
    GroundVeg::FT
    Tree::FT
    WallSun::FT
    WallShade::FT
    TotalGround::FT
    TotalCanyon::FT
    TotalUrban::FT
end

function DefaultRadiationFluxVariablesSubset(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, DefaultRadiationFluxVariablesSubset, Dict{String,Any}())
end

"""
    AlbedoOutput{FT<:AbstractFloat} <: AbstractModelVariables{FT}

Albedo values for urban areas.

# Fields
- `TotalUrban`: Albedo of the total urban area [-]
- `TotalCanyon`: Albedo of the total canyon area [-]
- `Roof`: Albedo of the total roof area [-]
"""
Base.@kwdef mutable struct AlbedoOutput{FT<:AbstractFloat} <: AbstractModelVariables{FT}
    TotalUrban::FT
    TotalCanyon::FT
    Roof::FT
end

function AlbedoOutput(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, AlbedoOutput, Dict{String,Any}())
end

"""
    RadiationFluxVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}

Container for all radiation variable components.

# Fields
- `SWRabs`: Absorbed shortwave radiation for different urban surfaces
- `SWRin`: Incoming shortwave radiation for different urban surfaces
- `SWRout`: Outgoing shortwave radiation for different urban surfaces
- `SWREB`: Shortwave radiation energy balance for different urban surfaces
- `LWRabs`: Absorbed longwave radiation for different urban surfaces
- `LWRin`: Incoming longwave radiation for different urban surfaces
- `LWRout`: Outgoing longwave radiation for different urban surfaces
- `LWREB`: Longwave radiation energy balance for different urban surfaces
- `AlbedoOutput`: Albedo values for urban areas
"""
Base.@kwdef struct RadiationFluxVariables{FT<:AbstractFloat} <: AbstractModelVariableSet{FT}
    SWRabs::AbsorbedRadiationFluxVariablesSubset{FT}
    SWRin::DefaultRadiationFluxVariablesSubset{FT}
    SWRout::DefaultRadiationFluxVariablesSubset{FT}
    SWREB::DefaultRadiationFluxVariablesSubset{FT}
    LWRabs::DefaultRadiationFluxVariablesSubset{FT}
    LWRin::DefaultRadiationFluxVariablesSubset{FT}
    LWRout::DefaultRadiationFluxVariablesSubset{FT}
    LWREB::DefaultRadiationFluxVariablesSubset{FT}
    AlbedoOutput::AlbedoOutput{FT}
end

function RadiationFluxVariables(::Type{FT}) where {FT<:AbstractFloat}
    return initialize(FT, RadiationFluxVariables, Dict{String,Any}())
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{RadiationFluxVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["SWRabs"] = AbsorbedRadiationFluxVariablesSubset(FT)
    processed["SWRin"] = DefaultRadiationFluxVariablesSubset(FT)
    processed["SWRout"] = DefaultRadiationFluxVariablesSubset(FT)
    processed["SWREB"] = DefaultRadiationFluxVariablesSubset(FT)
    processed["LWRabs"] = DefaultRadiationFluxVariablesSubset(FT)
    processed["LWRin"] = DefaultRadiationFluxVariablesSubset(FT)
    processed["LWRout"] = DefaultRadiationFluxVariablesSubset(FT)
    processed["LWREB"] = DefaultRadiationFluxVariablesSubset(FT)
    processed["AlbedoOutput"] = AlbedoOutput(FT)
    return processed
end

function ModelComponents.outputs_to_save(
    ::Type{RadiationFluxVariables}, ::Type{EssentialOutputs}
)
    return (:AlbedoOutput,)
end

function ModelComponents.outputs_to_save(
    ::Type{RadiationFluxVariables}, ::Type{ExtendedEnergyClimateOutputs}
)
    return (:SWRabs, :LWRabs)
end

function ModelComponents.outputs_to_save(
    ::Type{RadiationFluxVariables}, ::Type{ExtendedOutputs}
)
    return (:SWRin, :SWRout, :SWREB, :LWRin, :LWRout, :LWREB)
end
