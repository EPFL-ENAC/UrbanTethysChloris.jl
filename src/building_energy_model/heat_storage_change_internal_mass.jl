"""
    heat_storage_change_internal_mass(
        Tintmass::FT,
        Tintmasstm1::FT,
        ParThermalBuildFloor::NamedTuple,
        Geometry_m::NamedTuple,
        ParCalculation::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate heat storage change in building internal mass.

# Arguments
- `Tintmass`: Internal mass temperature [K]
- `Tintmasstm1`: Internal mass temperature at previous timestep [K]
- `ParThermalBuildFloor`: Thermal parameters for building floor
- `Geometry_m`: Building geometry parameters
- `ParCalculation`: Calculation parameters

# Returns
- `dS`: Energy storage change in internal mass [W/m²]
"""
function heat_storage_change_internal_mass(
    Tintmass::FT,
    Tintmasstm1::FT,
    ParThermalBuildFloor::NamedTuple,
    Geometry_m::NamedTuple,
    ParCalculation::NamedTuple,
) where {FT<:AbstractFloat}
    # Extract parameters
    Ts = Tintmass
    Ts_tm1 = Tintmasstm1
    dts = ParCalculation.dts
    cv_floor = ParThermalBuildFloor.cv_floor_IntMass
    cv_wall = ParThermalBuildFloor.cv_wall_IntMass
    dzFloor = ParThermalBuildFloor.dzFloor
    dzWall = ParThermalBuildFloor.dzWall
    FloorHeight = ParThermalBuildFloor.FloorHeight
    BuildingHeight = Geometry_m.Height_canyon
    BuildingWidth = Geometry_m.Width_roof

    # Calculate number of floor layers
    NrOfFloorLayers = max(round(Int, BuildingHeight/FloorHeight) - 1, 0)

    # Calculate internal mass volumes
    IntFloorVolume = NrOfFloorLayers * dzFloor * BuildingWidth
    IntWallVolume = dzWall * BuildingHeight

    # Normalize by building height
    IntTotal_dz = (IntFloorVolume + IntWallVolume) / BuildingHeight

    # Calculate average heat capacity
    cv_int_total = (
        cv_floor * IntFloorVolume / (IntFloorVolume + IntWallVolume) +
        cv_wall * IntWallVolume / (IntFloorVolume + IntWallVolume)
    )

    # Calculate heat storage change
    dS = cv_int_total * IntTotal_dz / dts * (Ts - Ts_tm1)

    return dS
end
