"""
    eb_solver_urban_climate_building_energy_model(
        TemperatureTot::Vector{FT},
        TempVec_ittm::NamedTuple,
        TempVecB_ittm::NamedTuple,
        Humidity_ittm::NamedTuple,
        MeteoData::NamedTuple,
        Int_ittm::NamedTuple,
        ExWater_ittm::NamedTuple,
        Vwater_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        SoilPotW_ittm::NamedTuple,
        CiCO2Leaf_ittm::NamedTuple,
        TempDamp_ittm::NamedTuple,
        ViewFactor::RayTracing.ViewFactor{FT},
        Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        WallLayers::NamedTuple,
        ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParInterceptionTree::NamedTuple,
        PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
        PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
        PropOpticalTree::ModelComponents.Parameters.SimpleOpticalProperties{FT},
        ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
        ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        SunPosition::NamedTuple,
        HumidityAtm::NamedTuple,
        Anthropogenic::NamedTuple,
        ParCalculation::NamedTuple,
        PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
        ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
        ParThermalBuildingFloor::ModelComponents.Parameters.ThermalBuilding{FT},
        ParWindows::ModelComponents.Parameters.WindowParameters{FT},
        BEM_on::Bool,
        RESPreCalc::Bool,
        fconvPreCalc::FT,
        fconv::FT,
        rsRoofPreCalc::NamedTuple,
        rsGroundPreCalc::NamedTuple,
        rsTreePreCalc::NamedTuple,
        HVACSchedule::NamedTuple
    ) where {FT<:AbstractFloat}

Solve energy balance for urban climate and building energy model.

# Arguments
- `TemperatureTot`: Temperature vector containing all surface temperatures
- `TempVec_ittm`: Temperature vectors at previous time step
- `TempVecB_ittm`: Building temperature vectors at previous time step
- `Humidity_ittm`: Humidity at previous time step
- Other parameters match those of eb_solver_roof, eb_solver_canyon and eb_solver_building

# Returns
- `Ytot`: Combined energy balance residuals from roof, canyon and building
"""
function eb_solver_urban_climate_building_energy_model(
    TemperatureTot::Vector{FT},
    TempVec_ittm::NamedTuple,
    TempVecB_ittm::NamedTuple,
    Humidity_ittm::NamedTuple,
    MeteoData::NamedTuple,
    Int_ittm::NamedTuple,
    ExWater_ittm::NamedTuple,
    Vwater_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    SoilPotW_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    TempDamp_ittm::NamedTuple,
    ViewFactor::RayTracing.ViewFactor{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    WallLayers::NamedTuple,
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParInterceptionTree::NamedTuple,
    PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    PropOpticalTree::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    ParThermalGround::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    ParThermalRoof::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SunPosition::NamedTuple,
    HumidityAtm::NamedTuple,
    Anthropogenic::NamedTuple,
    ParCalculation::NamedTuple,
    PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    ParThermalBuildingFloor::ModelComponents.Parameters.ThermalBuilding{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
    RESPreCalc::Bool,
    fconvPreCalc::FT,
    fconv::FT,
    rsRoofPreCalc::NamedTuple,
    rsGroundPreCalc::NamedTuple,
    rsTreePreCalc::NamedTuple,
    HVACSchedule::NamedTuple,
) where {FT<:AbstractFloat}

    # Extract temperature components
    TemperatureR = zeros(FT, 4)
    TemperatureR[1] = TemperatureTot[1]  # Troof_imp
    TemperatureR[2] = TemperatureTot[2]  # Troof_veg
    TemperatureR[3] = TemperatureTot[3]  # Troof_interior_imp
    TemperatureR[4] = TemperatureTot[4]  # Troof_interior_veg

    TemperatureC = zeros(FT, 10)
    TemperatureC[1] = TemperatureTot[5]   # Temperature ground impervious area
    TemperatureC[2] = TemperatureTot[6]   # Temperature ground bare area
    TemperatureC[3] = TemperatureTot[7]   # Temperature ground vegetated area
    TemperatureC[4] = TemperatureTot[8]   # Temperature sunlit area
    TemperatureC[5] = TemperatureTot[9]   # Temperature shaded area
    TemperatureC[6] = TemperatureTot[10]  # Temperature tree canopy
    TemperatureC[7] = TemperatureTot[11]  # Interior temperature sunlit wall
    TemperatureC[8] = TemperatureTot[12]  # Interior temperature shaded wall
    TemperatureC[9] = TemperatureTot[13]  # Temperature canyon
    TemperatureC[10] = TemperatureTot[14] # Specific humidity canyon

    TemperatureB = zeros(FT, 8)
    TemperatureB[1] = TemperatureTot[15]  # Temperature ceiling
    TemperatureB[2] = TemperatureTot[16]  # Temperature sunlit wall
    TemperatureB[3] = TemperatureTot[17]  # Temperature shaded wall
    TemperatureB[4] = TemperatureTot[18]  # Temperature windows
    TemperatureB[5] = TemperatureTot[19]  # Temperature ground
    TemperatureB[6] = TemperatureTot[20]  # Temperature building internal mass
    TemperatureB[7] = TemperatureTot[21]  # Temperature air
    TemperatureB[8] = TemperatureTot[22]  # Humidity air

    # Solve energy balance for roof surfaces
    Yroof, G2Roof = eb_solver_roof(
        TemperatureR,
        TemperatureB,
        TempVec_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        Gemeotry_m,
        FractionsRoof,
        ParSoilRoof,
        PropOpticalRoof,
        ParThermalRoof,
        ParVegRoof,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        BEM_on,
        RESPreCalc,
        rsRoofPreCalc,
    )

    # Solve energy balance for canyon surfaces
    Ycanyon, G2WallSun, G2WallShade, SWRabs_t, SWRinWsun, SWRinWshd = eb_solver_canyon(
        TemperatureC,
        TemperatureB,
        TempVec_ittm,
        Humidity_ittm,
        MeteoData,
        Int_ittm,
        ExWater_ittm,
        Vwater_ittm,
        Owater_ittm,
        SoilPotW_ittm,
        CiCO2Leaf_ittm,
        TempDamp_ittm,
        ViewFactor,
        Gemeotry_m,
        FractionsGround,
        WallLayers,
        ParSoilGround,
        ParInterceptionTree,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParThermalGround,
        ParThermalWall,
        ParVegGround,
        ParVegTree,
        SunPosition,
        HumidityAtm,
        Anthropogenic,
        ParCalculation,
        TempVecB_ittm,
        G2Roof,
        PropOpticalIndoors,
        ParHVAC,
        ParThermalBuildingFloor,
        ParWindows,
        BEM_on,
        RESPreCalc,
        fconvPreCalc,
        fconv,
        rsGroundPreCalc,
        rsTreePreCalc,
        HVACSchedule,
    )

    # Solve energy balance for building interior
    YBuildInt, _ = BuildingEnergyModel.eb_solver_building(
        TemperatureC,
        TemperatureB,
        TempVecB_ittm,
        TempVec_ittm,
        Humidity_ittm,
        MeteoData,
        SWRinWsun,
        SWRinWshd,
        G2Roof,
        G2WallSun,
        G2WallShade,
        TempDamp_ittm,
        SWRabs_t,
        Gemeotry_m,
        PropOpticalIndoors,
        ParHVAC,
        ParCalculation,
        ParThermalBuildingFloor,
        ParWindows,
        BEM_on,
        HVACSchedule,
    )

    # Combine all residuals
    Ytot = vcat(Yroof, Ycanyon, YBuildInt)

    return Ytot
end
