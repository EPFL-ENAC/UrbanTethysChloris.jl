"""
    f_solver_tot(
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
        TempVec_ittm2Ext::NamedTuple,
        Humidity_ittm2Ext::NamedTuple,
        TempVecB_ittm2Ext::NamedTuple,
        Meteo_ittm::NamedTuple,
        RESPreCalc::Bool,
        fconvPreCalc::FT,
        fconv::FT,
        rsRoofPreCalc::NamedTuple,
        rsGroundPreCalc::NamedTuple,
        rsTreePreCalc::NamedTuple,
        HVACSchedule::NamedTuple,
    ) where {FT<:AbstractFloat}

Nonlinear system solver for the coupled energy balance equations.

# Returns
- `T::Vector{FT}`: Solution vector containing temperatures
- `fval::Vector{FT}`: Function values at solution
- `exitflag::Int`: Exit flag indicating convergence status

# Temperature vector indices:
 1. Troof_imp: Impervious roof temperature
 2. Troof_veg: Vegetated roof temperature
 3. Troof_interior_imp: Interior impervious roof temperature
 4. Troof_interior_veg: Interior vegetated roof temperature
 5. TGroundImp: Ground impervious area temperature
 6. TGroundBare: Ground bare area temperature
 7. TGroundVeg: Ground vegetated area temperature
 8. TWallSun: Sunlit wall temperature
 9. TWallShade: Shaded wall temperature
10. TTree: Tree canopy temperature
11. TWallIntSun: Interior sunlit wall temperature
12. TWallIntShade: Interior shaded wall temperature
13. TCanyon: Canyon air temperature
14. qCanyon: Canyon specific humidity
15. TCeiling: Building ceiling temperature
16. TWallSunInt: Interior sunlit wall temperature
17. TWallShadeInt: Interior shaded wall temperature
18. TWindows: Window temperature
19. TGroundInt: Interior ground temperature
20. TIntMass: Internal mass temperature
21. TBin: Indoor air temperature
22. qBin: Indoor specific humidity
"""
function f_solver_tot(
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
    TempVec_ittm2Ext::NamedTuple,
    Humidity_ittm2Ext::NamedTuple,
    TempVecB_ittm2Ext::NamedTuple,
    Meteo_ittm::NamedTuple,
    RESPreCalc::Bool,
    fconvPreCalc::FT,
    fconv::FT,
    rsRoofPreCalc::NamedTuple,
    rsGroundPreCalc::NamedTuple,
    rsTreePreCalc::NamedTuple,
    HVACSchedule::NamedTuple;
    iterations::Int=500,
    f_tol::Real=1e-10,
) where {FT<:AbstractFloat}
    # Create lower bounds vector
    lb = fill(FT(243), 22)
    lb[14] = 0  # Lower bound for specific humidity
    lb[22] = 0  # Lower bound for specific humidity

    ub = fill(+Inf, 22)

    # Initialize storage for function values and exit flags
    # Initialize a vector of vector with typemax(FT) values

    fvals = fill(fill(typemax(FT), 22), 6)
    ran_simulation = [true, false, false, false, false, false]
    Ts = Vector{Vector{FT}}(undef, 6)
    exitflags = Vector{Bool}(undef, 6)

    # Attempt 1: Use previous timestep as initial guess
    if BEM_on
        TemperatureTot = [
            TempVec_ittm2Ext.TRoofImp[3],
            TempVec_ittm2Ext.TRoofVeg[3],
            TempVec_ittm2Ext.TRoofIntImp[3],
            TempVec_ittm2Ext.TRoofIntVeg[3],
            TempVec_ittm2Ext.TGroundImp[3],
            TempVec_ittm2Ext.TGroundBare[3],
            TempVec_ittm2Ext.TGroundVeg[3],
            TempVec_ittm2Ext.TWallSun[3],
            TempVec_ittm2Ext.TWallShade[3],
            TempVec_ittm2Ext.TTree[3],
            TempVec_ittm2Ext.TWallIntSun[3],
            TempVec_ittm2Ext.TWallIntShade[3],
            TempVec_ittm2Ext.TCanyon[3],
            Humidity_ittm2Ext.CanyonSpecific[3],
            TempVecB_ittm2Ext.Tceiling[3],
            TempVecB_ittm2Ext.Tinwallsun[3],
            TempVecB_ittm2Ext.Tinwallshd[3],
            TempVecB_ittm2Ext.Twindows[3],
            TempVecB_ittm2Ext.Tinground[3],
            TempVecB_ittm2Ext.Tintmass[3],
            TempVecB_ittm2Ext.Tbin[3],
            TempVecB_ittm2Ext.qbin[3],
        ]
    else
        TemperatureTot = [
            TempVec_ittm2Ext.TRoofImp[3],
            TempVec_ittm2Ext.TRoofVeg[3],
            TempVec_ittm2Ext.TRoofIntImp[3],
            TempVec_ittm2Ext.TRoofIntVeg[3],
            TempVec_ittm2Ext.TGroundImp[3],
            TempVec_ittm2Ext.TGroundBare[3],
            TempVec_ittm2Ext.TGroundVeg[3],
            TempVec_ittm2Ext.TWallSun[3],
            TempVec_ittm2Ext.TWallShade[3],
            TempVec_ittm2Ext.TTree[3],
            TempVec_ittm2Ext.TWallIntSun[3],
            TempVec_ittm2Ext.TWallIntShade[3],
            TempVec_ittm2Ext.TCanyon[3],
            Humidity_ittm2Ext.CanyonSpecific[3],
            FT(273.15),
            FT(273.15),
            FT(273.15),
            FT(273.15),
            FT(273.15),
            FT(273.15),
            FT(273.15),
            FT(0),
        ]
    end

    # TODO: set 273.15 as global constant

    # Define target function for optimization
    function target_fun(TemperatureTot)
        return eb_solver_urban_climate_building_energy_model(
            TemperatureTot,
            TempVec_ittm,
            TempVecB_ittm,
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
            FractionsRoof,
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
            ParSoilRoof,
            PropOpticalRoof,
            ParThermalRoof,
            ParVegRoof,
            SunPosition,
            HumidityAtm,
            Anthropogenic,
            ParCalculation,
            PropOpticalIndoors,
            ParHVAC,
            ParThermalBuildingFloor,
            ParWindows,
            BEM_on,
            RESPreCalc,
            fconvPreCalc,
            fconv,
            rsRoofPreCalc,
            rsGroundPreCalc,
            rsTreePreCalc,
            HVACSchedule,
        )
    end

    function lsqnonlin(x0::Vector{FT})
        results = optimize(
            target_fun,
            x0,
            LevenbergMarquardt();
            lower=lb,
            upper=ub,
            iterations=iterations,
            f_tol=f_tol,
        )

        T = results.minimizer
        fval = target_fun(T)
        exitflag = results.converged

        return T, fval, exitflag
    end

    Ts[1], fvals[1], exitflags[1] = lsqnonlin(TemperatureTot)

    if sum(abs.(fvals[1])) > 0.1
        return Ts[1], fvals[1], exitflags[1]
    end

    # If first attempt failed, try with different initial conditions
    if maximum(abs.(fvals[1])) > 0.1
        # Additional attempts with different initial conditions
        # Attempt 2: Perturbed initial conditions based on solar radiation
        if Meteo_ittm.SWRin[2] > Meteo_ittm.SWRin[1] && Meteo_ittm.Rain[2] < 1
            r = 3 .* rand(FT, 22)
            r[13] = rand(FT)
            r[21] = rand(FT)
            r[14] = -1e-4
            r[22] = -1e-4
        elseif Meteo_ittm.SWRin[2] <= Meteo_ittm.SWRin[1]
            r = -3 .* rand(FT, 22)
            r[13] = -rand(FT)
            r[21] = -rand(FT)
            r[14] = -1e-4
            r[22] = -1e-4
        else
            r = (-1 + 2 * rand(FT, 22)) .* 3
            r[14] = -1e-4
            r[22] = -1e-4
        end

        if BEM_on
            TemperatureTot = [
                TempVec_ittm2Ext.TRoofImp[3] + r[1],
                TempVec_ittm2Ext.TRoofVeg[3] + r[2],
                TempVec_ittm2Ext.TRoofIntImp[3] + r[3],
                TempVec_ittm2Ext.TRoofIntVeg[3] + r[4],
                TempVec_ittm2Ext.TGroundImp[3] + r[5],
                TempVec_ittm2Ext.TGroundBare[3] + r[6],
                TempVec_ittm2Ext.TGroundVeg[3] + r[7],
                TempVec_ittm2Ext.TWallSun[3] + r[8],
                TempVec_ittm2Ext.TWallShade[3] + r[9],
                TempVec_ittm2Ext.TTree[3] + r[10],
                TempVec_ittm2Ext.TWallIntSun[3] + r[11],
                TempVec_ittm2Ext.TWallIntShade[3] + r[12],
                TempVec_ittm2Ext.TCanyon[3] + r[13],
                min(Humidity_ittm2Ext.CanyonSpecific[3] + r[14], FT(0)),
                TempVecB_ittm2Ext.Tceiling[3] + r[15],
                TempVecB_ittm2Ext.Tinwallsun[3] + r[16],
                TempVecB_ittm2Ext.Tinwallshd[3] + r[17],
                TempVecB_ittm2Ext.Twindows[3] + r[18],
                TempVecB_ittm2Ext.Tinground[3] + r[19],
                TempVecB_ittm2Ext.Tintmass[3] + r[20],
                TempVecB_ittm2Ext.Tbin[3] + r[21],
                min(FT(0), Humidity_ittm2Ext.qbin[3] + r[22]),
            ]
        else
            TemperatureTot = [
                TempVec_ittm2Ext.TRoofImp[3] + r[1],
                TempVec_ittm2Ext.TRoofVeg[3] + r[2],
                TempVec_ittm2Ext.TRoofIntImp[3] + r[3],
                TempVec_ittm2Ext.TRoofIntVeg[3] + r[4],
                TempVec_ittm2Ext.TGroundImp[3] + r[5],
                TempVec_ittm2Ext.TGroundBare[3] + r[6],
                TempVec_ittm2Ext.TGroundVeg[3] + r[7],
                TempVec_ittm2Ext.TWallSun[3] + r[8],
                TempVec_ittm2Ext.TWallShade[3] + r[9],
                TempVec_ittm2Ext.TTree[3] + r[10],
                TempVec_ittm2Ext.TWallIntSun[3] + r[11],
                TempVec_ittm2Ext.TWallIntShade[3] + r[12],
                TempVec_ittm2Ext.TCanyon[3] + r[13],
                Humidity_ittm2Ext.CanyonSpecific[3] + r[14],
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                0.0,
            ]
        end

        Ts[2], fvals[2], exitflags[2] = lsqnonlin(TemperatureTot)

        ran_simulation[2] = true
    end

    if ran_simulation[2] && maximum(abs.(fvals[2])) > 0.1
        TT = MeteoData.Tatm
        # Attempt 3: Retry with atmospheric forcing
        if BEM_on
            TemperatureTot = [
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                MeteoData.q_atm,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                MeteoData.q_atm,
            ]
        else
            TemperatureTot = [
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                TT,
                MeteoData.q_atm,
                273.15,
                273.15,
                273.15,
                273.15,
                273.15,
                273.15,
                273.15,
                0,
            ];
        end

        Ts[3], fvals[3], exitflags[3] = lsqnonlin(TemperatureTot)

        ran_simulation[2] = true
    end

    if ran_simulation[3] && maximum(abs.(fvals[3])) > 0.01
        # Attempt 4: Retry with two steps extrapolated
        if BEM_on
            TemperatureTot = [
                TempVec_ittm2Ext.TRoofImp[4],
                TempVec_ittm2Ext.TRoofVeg[4],
                TempVec_ittm2Ext.TRoofIntImp[4],
                TempVec_ittm2Ext.TRoofIntVeg[4],
                TempVec_ittm2Ext.TGroundImp[4],
                TempVec_ittm2Ext.TGroundBare[4],
                TempVec_ittm2Ext.TGroundVeg[4],
                TempVec_ittm2Ext.TWallSun[4],
                TempVec_ittm2Ext.TWallShade[4],
                TempVec_ittm2Ext.TTree[4],
                TempVec_ittm2Ext.TWallIntSun[4],
                TempVec_ittm2Ext.TWallIntShade[4],
                TempVec_ittm2Ext.TCanyon[4],
                Humidity_ittm2Ext.CanyonSpecific[4],
                TempVecB_ittm2Ext.Tceiling[4],
                TempVecB_ittm2Ext.Tinwallsun[4],
                TempVecB_ittm2Ext.Tinwallshd[4],
                TempVecB_ittm2Ext.Twindows[4],
                TempVecB_ittm2Ext.Tinground[4],
                TempVecB_ittm2Ext.Tintmass[4],
                TempVecB_ittm2Ext.Tbin[4],
                TempVecB_ittm2Ext.qbin[4],
            ]
        else
            TemperatureTot = [
                TempVec_ittm2Ext.TRoofImp[4],
                TempVec_ittm2Ext.TRoofVeg[4],
                TempVec_ittm2Ext.TRoofIntImp[4],
                TempVec_ittm2Ext.TRoofIntVeg[4],
                TempVec_ittm2Ext.TGroundImp[4],
                TempVec_ittm2Ext.TGroundBare[4],
                TempVec_ittm2Ext.TGroundVeg[4],
                TempVec_ittm2Ext.TWallSun[4],
                TempVec_ittm2Ext.TWallShade[4],
                TempVec_ittm2Ext.TTree[4],
                TempVec_ittm2Ext.TWallIntSun[4],
                TempVec_ittm2Ext.TWallIntShade[4],
                TempVec_ittm2Ext.TCanyon[4],
                Humidity_ittm2Ext.CanyonSpecific[4],
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(0),
            ]
        end

        Ts[4], fvals[4], exitflags[4] = lsqnonlin(TemperatureTot)

        ran_simulation[4] = true
    end

    if ran_simulation[4] && maximum(abs.(fvals[4])) > 0.01
        # Attempt 5: Retry with previous time step

        if BEM_on
            TemperatureTot = [
                TempVec_ittm2Ext.TRoofImp[2],
                TempVec_ittm2Ext.TRoofVeg[2],
                TempVec_ittm2Ext.TRoofIntImp[2],
                TempVec_ittm2Ext.TRoofIntVeg[2],
                TempVec_ittm2Ext.TGroundImp[2],
                TempVec_ittm2Ext.TGroundBare[2],
                TempVec_ittm2Ext.TGroundVeg[2],
                TempVec_ittm2Ext.TWallSun[2],
                TempVec_ittm2Ext.TWallShade[2],
                TempVec_ittm2Ext.TTree[2],
                TempVec_ittm2Ext.TWallIntSun[2],
                TempVec_ittm2Ext.TWallIntShade[2],
                TempVec_ittm2Ext.TCanyon[2],
                Humidity_ittm2Ext.CanyonSpecific[2],
                TempVecB_ittm2Ext.Tceiling[2],
                TempVecB_ittm2Ext.Tinwallsun[2],
                TempVecB_ittm2Ext.Tinwallshd[2],
                TempVecB_ittm2Ext.Twindows[2],
                TempVecB_ittm2Ext.Tinground[2],
                TempVecB_ittm2Ext.Tintmass[2],
                TempVecB_ittm2Ext.Tbin[2],
                TempVecB_ittm2Ext.qbin[2],
            ]
        else
            TemperatureTot = [
                TempVec_ittm2Ext.TRoofImp[2],
                TempVec_ittm2Ext.TRoofVeg[2],
                TempVec_ittm2Ext.TRoofIntImp[2],
                TempVec_ittm2Ext.TRoofIntVeg[2],
                TempVec_ittm2Ext.TGroundImp[2],
                TempVec_ittm2Ext.TGroundBare[2],
                TempVec_ittm2Ext.TGroundVeg[2],
                TempVec_ittm2Ext.TWallSun[2],
                TempVec_ittm2Ext.TWallShade[2],
                TempVec_ittm2Ext.TTree[2],
                TempVec_ittm2Ext.TWallIntSun[2],
                TempVec_ittm2Ext.TWallIntShade[2],
                TempVec_ittm2Ext.TCanyon[2],
                Humidity_ittm2Ext.CanyonSpecific[2],
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(0),
            ]
        end

        Ts[5], fvals[5], exitflags[5] = lsqnonlin(TemperatureTot)

        ran_simulation[5] = true
    end

    if ran_simulation[5] && maximum(abs.(fvals[5])) > 0.1
        # Attempt 6: Retry randomly

        TT = MeteoData.Tatm

        r = (-1 + 2 * rand(FT, 22)) .* 3

        if BEM_on
            TemperatureTot = [
                TT + r[1] * (TempVec_ittm.TRoofImp - 273.15),
                TT + r[2] * (TempVec_ittm.TRoofVeg - 273.15),
                TT + r[3] * (TempVec_ittm.TRoofIntImp - 273.15),
                TT + r[4] * (TempVec_ittm.TRoofIntVeg - 273.15),
                TT + r[5] * (TempVec_ittm.TGroundImp - 273.15),
                TT + r[6] * (TempVec_ittm.TGroundBare - 273.15),
                TT + r[7] * (TempVec_ittm.TGroundVeg - 273.15),
                TT + r[8] * (TempVec_ittm.TWallSun - 273.15),
                TT + r[9] * (TempVec_ittm.TWallShade - 273.15),
                TT + r[10] * (TempVec_ittm.TTree - 273.15),
                TT + r[11] * (TempVec_ittm.TWallIntSun - 273.15),
                TT + r[12] * (TempVec_ittm.TWallIntShade - 273.15),
                TT + r[13] * (TempVec_ittm.TCanyon - 273.15),
                MeteoData.q_atm,
                TT + r[15] * (TempVec_ittm.Tceiling - 273.15),
                TT + r[16] * (TempVec_ittm.Tinwallsun - 273.15),
                TT + r[17] * (TempVec_ittm.Tinwallshd - 273.15),
                TT + r[18] * (TempVec_ittm.Twindows - 273.15),
                TT + r[19] * (TempVec_ittm.Tinground - 273.15),
                TT + r[20] * (TempVec_ittm.Tintmass - 273.15),
                TT + r[21] * (TempVec_ittm.Tbin - 273.15),
                MeteoData.q_atm,
            ]
        else
            TemperatureTot = [
                TT + r[1] * (TempVec_ittm.TRoofImp - 273.15),
                TT + r[2] * (TempVec_ittm.TRoofVeg - 273.15),
                TT + r[3] * (TempVec_ittm.TRoofIntImp - 273.15),
                TT + r[4] * (TempVec_ittm.TRoofIntVeg - 273.15),
                TT + r[5] * (TempVec_ittm.TGroundImp - 273.15),
                TT + r[6] * (TempVec_ittm.TGroundBare - 273.15),
                TT + r[7] * (TempVec_ittm.TGroundVeg - 273.15),
                TT + r[8] * (TempVec_ittm.TWallSun - 273.15),
                TT + r[9] * (TempVec_ittm.TWallShade - 273.15),
                TT + r[10] * (TempVec_ittm.TTree - 273.15),
                TT + r[11] * (TempVec_ittm.TWallIntSun - 273.15),
                TT + r[12] * (TempVec_ittm.TWallIntShade - 273.15),
                TT + r[13] * (TempVec_ittm.TCanyon - 273.15),
                MeteoData.q_atm,
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                FT(273.15),
                0.0,
            ]
        end

        Ts[6], fvals[6], exitflags[6] = lsqnonlin(TemperatureTot)

        ran_simulation[6] = true
    end

    # Select best solution

    best_idx = findmin(x -> maximum(abs.(x)), fvals)[2]
    T = Ts[best_idx]
    fval = fvals[best_idx]
    exitflag = exitflags[best_idx]

    return T, fval, exitflag
end
