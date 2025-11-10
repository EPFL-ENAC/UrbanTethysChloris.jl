"""
    AbsorbedRadiationFluxVariablesSubset{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct AbsorbedRadiationFluxVariablesSubset{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    SWRabsRoofImp::Array{FT,N}
    RoofVeg::Array{FT,N}
    TotalRoof::Array{FT,N}
    GroundImp::Array{FT,N}
    GroundBare::Array{FT,N}
    GroundVeg::Array{FT,N}
    Tree::Array{FT,N}
    WallSun::Array{FT,N}
    WallShade::Array{FT,N}
    TotalGround::Array{FT,N}
    TotalCanyon::Array{FT,N}
    TotalUrban::Array{FT,N}
    WallSunExt::Array{FT,N}
    WallShadeExt::Array{FT,N}
    WallSunTransmitted::Array{FT,N}
    WallShadeTransmitted::Array{FT,N}
end

function initialize_absorbed_radiation_flux_variables(
    ::Type{FT}, ::TimeSlice
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        AbsorbedRadiationFluxVariablesSubset,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSlice())),
    )
end

function initialize_absorbed_radiation_flux_variables(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        AbsorbedRadiationFluxVariablesSubset,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
    )
end

"""
    DefaultRadiationFluxVariablesSubset{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

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
Base.@kwdef mutable struct DefaultRadiationFluxVariablesSubset{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    RoofImp::Array{FT,N}
    RoofVeg::Array{FT,N}
    TotalRoof::Array{FT,N}
    GroundImp::Array{FT,N}
    GroundBare::Array{FT,N}
    GroundVeg::Array{FT,N}
    Tree::Array{FT,N}
    WallSun::Array{FT,N}
    WallShade::Array{FT,N}
    TotalGround::Array{FT,N}
    TotalCanyon::Array{FT,N}
    TotalUrban::Array{FT,N}
end

function initialize_default_radiation_flux_variables(
    ::Type{FT}, ::TimeSlice
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        DefaultRadiationFluxVariablesSubset,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSlice())),
    )
end

function initialize_default_radiation_flux_variables(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        DefaultRadiationFluxVariablesSubset,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
    )
end

"""
    AlbedoOutput{FT<:AbstractFloat, N} <: Abstract1PModelVariables{FT,N}

Albedo values for urban areas.

# Fields
- `TotalUrban`: Albedo of the total urban area [-]
- `TotalCanyon`: Albedo of the total canyon area [-]
- `Roof`: Albedo of the total roof area [-]
"""
Base.@kwdef mutable struct AlbedoOutput{FT<:AbstractFloat,N} <:
                           Abstract1PModelVariables{FT,N}
    TotalUrban::Array{FT,N}
    TotalCanyon::Array{FT,N}
    Roof::Array{FT,N}
end

function initialize_albedo_output(::Type{FT}, ::TimeSlice) where {FT<:AbstractFloat}
    return initialize(
        FT, AlbedoOutput, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function initialize_albedo_output(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT, AlbedoOutput, Dict{String,Any}(), (FT, dimension_value(TimeSeries())), hours
    )
end

"""
    RadiationFluxVariables{FT<:AbstractFloat, N} <: Abstract1PModelVariablesSet{FT, N}

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
Base.@kwdef struct RadiationFluxVariables{FT<:AbstractFloat,N} <:
                   Abstract1PModelVariablesSet{FT,N}
    SWRabs::AbsorbedRadiationFluxVariablesSubset{FT,N}
    SWRin::DefaultRadiationFluxVariablesSubset{FT,N}
    SWRout::DefaultRadiationFluxVariablesSubset{FT,N}
    SWREB::DefaultRadiationFluxVariablesSubset{FT,N}
    LWRabs::DefaultRadiationFluxVariablesSubset{FT,N}
    LWRin::DefaultRadiationFluxVariablesSubset{FT,N}
    LWRout::DefaultRadiationFluxVariablesSubset{FT,N}
    LWREB::DefaultRadiationFluxVariablesSubset{FT,N}
    AlbedoOutput::AlbedoOutput{FT,N}
end

function initialize_radiation_flux_variables(
    ::Type{FT}, ::TimeSlice
) where {FT<:AbstractFloat}
    return initialize(
        FT, RadiationFluxVariables, Dict{String,Any}(), (FT, dimension_value(TimeSlice()))
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT}, ::Type{RadiationFluxVariables}, data::Dict{String,Any}, params::Tuple
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["SWRabs"] = initialize_absorbed_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["SWRin"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["SWRout"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["SWREB"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["LWRabs"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["LWRin"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["LWRout"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["LWREB"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2])
    )
    processed["AlbedoOutput"] = initialize_albedo_output(FT, dimensionality_type(params[2]))
    return processed
end

function initialize_radiation_flux_variables(
    ::Type{FT}, ::TimeSeries, hours::Int
) where {FT<:AbstractFloat}
    return initialize(
        FT,
        RadiationFluxVariables,
        Dict{String,Any}(),
        (FT, dimension_value(TimeSeries())),
        hours,
    )
end

function TethysChlorisCore.preprocess_fields(
    ::Type{FT},
    ::Type{RadiationFluxVariables},
    data::Dict{String,Any},
    params::Tuple,
    hours::Int,
) where {FT<:AbstractFloat}
    processed = Dict{String,Any}()
    processed["SWRabs"] = initialize_absorbed_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["SWRin"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["SWRout"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["SWREB"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["LWRabs"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["LWRin"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["LWRout"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["LWREB"] = initialize_default_radiation_flux_variables(
        FT, dimensionality_type(params[2]), hours
    )
    processed["AlbedoOutput"] = initialize_albedo_output(
        FT, dimensionality_type(params[2]), hours
    )
    return processed
end
