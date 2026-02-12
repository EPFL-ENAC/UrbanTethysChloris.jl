function run_simulation(
    model::Model{FT},
    forcing::UrbanTethysChloris.ModelComponents.ForcingInputSet{FT,1};
    NN::Signed=nothing,
    mc_sample_size::Int=1000,
    n_rays::Int=200,
    RESPreCalc::Bool=true,
    fconvPreCalc::Bool=true,
    BEM_on::Bool=true,
    WallLayers::NamedTuple=(; dz1_wall=FT(0.1), dz2_wall=FT(0.1)),
    ParInterceptionTree::NamedTuple=(; Sp_In=FT(0.2)),
    ViewFactors::Union{
        Tuple{RayTracing.ViewFactor{FT},RayTracing.ViewFactorPoint{FT}},Nothing
    }=nothing,
    O33::NamedTuple=(roof=FT(0.0), ground=FT(0.0)),
    output_level::UrbanTethysChloris.ModelComponents.AbstractOutputsToSave=plot_outputs,
) where {FT<:AbstractFloat}

    # view factor
    if isnothing(ViewFactors)
        ViewFactor, ViewFactorPoint = RayTracing.view_factors_canyon(
            model.parameters.urbangeometry, model.parameters.person, mc_sample_size, n_rays
        )
    else
        ViewFactor, ViewFactorPoint = ViewFactors
    end

    # initialize Vwater and Owater
    # TODO: initialize in initialize! function?
    if model.parameters.surfacefractions.roof.fimp == 1
        model.variables.waterflux.Vwater.VRoofSoilVeg[:] .= 0
        model.variables.waterflux.Owater.OwRoofSoilVeg[:] .= 0
    end

    if model.parameters.surfacefractions.ground.fimp == 1
        model.variables.waterflux.Vwater.VGroundSoilImp[:] .= 0
        model.variables.waterflux.Owater.OGroundSoilImp[:] .= 0
        model.variables.waterflux.Vwater.VGroundSoilBare[:] .= 0
        model.variables.waterflux.Owater.OGroundSoilBare[:] .= 0
        model.variables.waterflux.Vwater.VGroundSoilVeg[:] .= 0
        model.variables.waterflux.Owater.OGroundSoilVeg[:] .= 0
    end

    # store OwaterInitial for back-computation
    OwaterInitial = deepcopy(model.variables.waterflux.Owater)

    # Hard-code ittm2Ext with values from first iteration, keep as named tuple for dev

    # Pre-calculate for faster numerical solution

    # Meteo_ittm
    ParCalculation = (;
        dth=1.0,
        dts=3600,
        row=1000.0,
        cp_atm=model.forcing.meteorological.cp_atm,
        rho_atm=model.forcing.meteorological.rho_atm,
    )

    Humidity_t = deepcopy(model.variables.humidity.Humidity)

    Ttot = nothing
    Yroof, Ycanyon, YBuildInt = nothing, nothing, nothing

    # TODO: create a modeltm1, modeltm2 structures, to stop all the _ittm and _ittm2Ext
    TempVec_ittm = deepcopy(model.variables.temperature.tempvec)
    Humidity_ittm = deepcopy(model.variables.humidity.Humidity)
    TempVecB_ittm = deepcopy(model.variables.buildingenergymodel.TempVecB)
    Int_ittm = deepcopy(model.variables.waterflux.Interception)
    ExWater_ittm = deepcopy(model.variables.waterflux.ExWater)
    Vwater_ittm = deepcopy(model.variables.waterflux.Vwater)
    Owater_ittm = deepcopy(model.variables.waterflux.Owater)
    SoilPotW_ittm = deepcopy(model.variables.waterflux.SoilPotW)
    CiCO2Leaf_ittm = deepcopy(model.variables.waterflux.CiCO2Leaf)
    Runon_ittm = deepcopy(model.variables.waterflux.Runon)
    TempDamp_ittm = deepcopy(model.variables.temperature.tempdamp)
    Qinlat_ittm = deepcopy(model.variables.waterflux.Qinlat)
    RES_ittm = deepcopy(model.variables.environmentalconditions.resistance)

    TempVec_ittm2Ext = ExtrapolatedTempVec(model.variables.temperature.tempvec)
    Humidity_ittm2Ext = ExtrapolatedHumidity(model.variables.humidity.Humidity)
    TempVecB_ittm2Ext = ExtrapolatedTempVecB(model.variables.buildingenergymodel.TempVecB)
    Meteo_ittm = Meteotm1(model.forcing.meteorological)

    SWRout_t = Radiation.RadiationFluxes(FT)
    LWRout_t = Radiation.RadiationFluxes(FT)
    T2m = FT(NaN)
    RH_T2m = FT(NaN)

    results_dict, accessors = prepare_results(typeof(output_level), model, NN)

    results_dict[:OwaterInitial] = Dict{Symbol,Array}(
        :OwRoofSoilVeg => OwaterInitial.OwRoofSoilVeg,
        :OwGroundSoilImp => OwaterInitial.OwGroundSoilImp,
        :OwGroundSoilBare => OwaterInitial.OwGroundSoilBare,
        :OwGroundSoilVeg => OwaterInitial.OwGroundSoilVeg,
        :OwGroundSoilTot => OwaterInitial.OwGroundSoilTot,
    )

    for i in 1:NN
        @info "Starting iteration $i / $NN"

        if i > 1
            extrapolate!(TempVec_ittm2Ext, model.variables.temperature.tempvec, i)
            extrapolate!(Humidity_ittm2Ext, model.variables.humidity.Humidity, i)
            extrapolate!(TempVecB_ittm2Ext, model.variables.buildingenergymodel.TempVecB, i)
            # TODO: bring HumidityAtm back to the forcing inputs
            # TODO rename as "apply forcing"
            update!(model.variables.humidity.Humidity, model.forcing.meteorological)
            update!(model.variables.temperature.tempvec, model.forcing.meteorological)
            update!(Meteo_ittm, model.forcing.meteorological)
        end

        if RESPreCalc || fconvPreCalc
            fconv, rsRoofPreCalc, rsGroundPreCalc, rsTreePreCalc = Resistance.precalculate_for_faster_numerical_solution(
                model,
                TempVec_ittm,
                Humidity_ittm,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                RES_ittm,
                i,
                1,
                ViewFactor,
                BEM_on,
            )
        else
            fconv = FT(NaN)
            rsRoofPreCalc = (;)
            rsGroundPreCalc = (;)
            rsTreePreCalc = (;)
        end

        ParHVAC, ParHVACorig = BuildingEnergyModel.ac_heating_turn_on_off(
            model, TempVecB_ittm, TempVec_ittm, Humidity_ittm, BEM_on
        )

        EnergyUse = (;);

        for HVACittm in 1:2
            if BEM_on && HVACittm == 2
                if !ParHVACorig.ACon && !ParHVACorig.Heatingon
                    continue
                end

                if EnergyUse.EnergyForAC_H>-1e-6 &&
                    EnergyUse.EnergyForAC_LE>-1e-6 &&
                    EnergyUse.EnergyForHeating>-1e-6
                    if ParHVAC.ACon &&
                        round(
                            model.variables.buildingenergymodel.TempVecB.Tbin; digits=4
                        )<(ParHVAC.TsetpointCooling+0.01) &&
                        round(
                            model.variables.buildingenergymodel.TempVecB.qbin; digits=8
                        )<(ParHVAC.q_RHspCooling+1e-6)
                        continue
                    elseif ParHVAC.Heatingon &&
                        round(
                        model.variables.buildingenergymodel.TempVecB.Tbin; digits=4
                    )>(ParHVAC.TsetpointHeating-0.01)
                        continue
                    end
                end

                ParHVAC = update_hvac_parameters(
                    ParHVACorig,
                    ParHVAC,
                    EnergyUse,
                    model.variables.buildingenergymodel.TempVecB.Tbin,
                    model.variables.buildingenergymodel.TempVecB.qbin,
                )
            end
            Ttot = f_solver_tot!(
                model,
                TempVec_ittm,
                TempVecB_ittm,
                Humidity_ittm,
                Int_ittm,
                ExWater_ittm,
                Vwater_ittm,
                Owater_ittm,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                TempDamp_ittm,
                ViewFactor,
                WallLayers,
                ParInterceptionTree,
                ParCalculation,
                ParHVAC,
                BEM_on,
                TempVec_ittm2Ext,
                Humidity_ittm2Ext,
                TempVecB_ittm2Ext,
                Meteo_ittm,
                RESPreCalc,
                fconvPreCalc,
                fconv,
                rsRoofPreCalc,
                rsGroundPreCalc,
                rsTreePreCalc,
            )

            update!(model.variables.temperature.tempvec, Ttot)
            update!(model.variables.humidity.Humidity, Ttot)
            update!(model.variables.buildingenergymodel.TempVecB, Ttot)

            TR = roof_temperature(model.variables.temperature.tempvec)
            TC = canyon_temperature(
                model.variables.temperature.tempvec, model.variables.humidity.Humidity
            )
            TB = building_temperature(model.variables.buildingenergymodel.TempVecB)

            # TODO: remove all EB, WB and Yroof from the list of outputs, they are already modified in-place
            G2Roof, Yroof = eb_wb_roof!(
                model,
                TR,
                TB,
                TempVec_ittm,
                Int_ittm,
                ExWater_ittm,
                Vwater_ittm,
                Owater_ittm,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                Runon_ittm,
                ParCalculation,
                BEM_on,
                RESPreCalc,
                rsRoofPreCalc,
            )

            SWRout_t, SWRabs_t, LWRout_t, G2WallSun, G2WallShade, Ycanyon, T2m, RH_T2m = eb_wb_canyon!(
                model,
                TC,
                TB,
                TempVec_ittm,
                Humidity_ittm,
                TempVecB_ittm,
                Int_ittm,
                ExWater_ittm,
                Vwater_ittm,
                Owater_ittm,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                TempDamp_ittm,
                Runon_ittm,
                Qinlat_ittm,
                ViewFactor,
                WallLayers,
                ParInterceptionTree,
                ParCalculation,
                G2Roof,
                ParHVAC,
                BEM_on,
                RESPreCalc,
                fconvPreCalc,
                fconv,
                rsGroundPreCalc,
                rsTreePreCalc,
            )

            SWRinWsun = SWRabs_t.WallSun
            SWRinWshd = SWRabs_t.WallShade

            EnergyUse, YBuildInt = BuildingEnergyModel.eb_solver_building_output!(
                model,
                TC,
                TB,
                TempVecB_ittm,
                TempVec_ittm,
                Humidity_ittm,
                TempDamp_ittm,
                SWRinWsun,
                SWRinWshd,
                G2Roof,
                G2WallSun,
                G2WallShade,
                SWRabs_t,
                ParHVAC,
                ParCalculation,
                BEM_on,
            )
        end

        MeanRadiantTemperature.mean_radiant_temperature!(
            model, SWRout_t, LWRout_t, ViewFactorPoint
        )
        # TODO: check whether we should be using the hour as a float (e.g. 10.5 for 10:30) or
        # the hour as an integer (10 for 10:30).

        Resistance.wind_profile_point_output!(model)

        OutdoorThermalComfort.utci_approx!(model)

        # Assign outputs
        model.variables.energybalance.Solver.YfunctionOutput = vcat(
            Yroof, Ycanyon, YBuildInt
        )

        # TODO: implement energy balance check script as function

        # tempvec
        model.variables.temperature.tempvec.T2m = T2m

        ModelComponents.fix_soil_moisture!(
            model.variables.waterflux.Owater,
            model.parameters.soil.roof,
            model.parameters.soil.ground,
            O33,
        )

        update!(TempVecB_ittm, model.variables.buildingenergymodel.TempVecB)
        update!(TempVec_ittm, model.variables.temperature.tempvec)
        update!(Humidity_ittm, model.variables.humidity.Humidity)
        update!(Int_ittm, model.variables.waterflux.Interception)
        update!(ExWater_ittm, model.variables.waterflux.ExWater)
        update!(Vwater_ittm, model.variables.waterflux.Vwater)
        update!(Owater_ittm, model.variables.waterflux.Owater)
        update!(SoilPotW_ittm, model.variables.waterflux.SoilPotW)
        update!(CiCO2Leaf_ittm, model.variables.waterflux.CiCO2Leaf)
        update!(TempDamp_ittm, model.variables.temperature.tempdamp)
        update!(Runon_ittm, model.variables.waterflux.Runon)
        update!(Qinlat_ittm, model.variables.waterflux.Qinlat)
        update!(RES_ittm, model.variables.environmentalconditions.resistance)

        urban_averages!(model)

        assign_results!(results_dict, accessors, model, i)

        # Update forcing parameters for the next step
        model.forcing = forcing[i + 1]
    end

    return results_dict, ViewFactor, ViewFactorPoint
end

function roof_temperature(
    tempvec::ModelComponents.ModelVariables.TempVec{FT}
) where {FT<:AbstractFloat}
    TR = zeros(FT, 4)

    TR[1] = tempvec.TRoofImp
    TR[2] = tempvec.TRoofVeg
    TR[3] = tempvec.TRoofIntImp
    TR[4] = tempvec.TRoofIntVeg

    return TR
end

function canyon_temperature(
    tempvec::ModelComponents.ModelVariables.TempVec{FT},
    humidity::ModelComponents.ModelVariables.Humidity{FT},
) where {FT<:AbstractFloat}
    TC = zeros(FT, 10)

    TC[1] = tempvec.TGroundImp
    TC[2] = tempvec.TGroundBare
    TC[3] = tempvec.TGroundVeg
    TC[4] = tempvec.TWallSun
    TC[5] = tempvec.TWallShade
    TC[6] = tempvec.TTree
    TC[7] = tempvec.TWallIntSun
    TC[8] = tempvec.TWallIntShade
    TC[9] = tempvec.TCanyon
    TC[10] = humidity.CanyonSpecific

    return TC
end

function building_temperature(
    tempvecB::ModelComponents.ModelVariables.TempVecB{FT}
) where {FT<:AbstractFloat}
    TB = zeros(FT, 8)

    TB[1] = tempvecB.Tceiling
    TB[2] = tempvecB.Tinwallsun
    TB[3] = tempvecB.Tinwallshd
    TB[4] = tempvecB.Twindows
    TB[5] = tempvecB.Tinground
    TB[6] = tempvecB.Tintmass
    TB[7] = tempvecB.Tbin
    TB[8] = tempvecB.qbin

    return TB
end

function urban_average(
    roof::FT,
    canyon::FT,
    urbangeometry::ModelComponents.Parameters.UrbanGeometryParameters{FT},
) where {FT<:AbstractFloat}
    froof = urbangeometry.wroof_norm
    fcanyon = urbangeometry.wcanyon_norm

    return roof * froof + canyon * fcanyon
end

function urban_averages!(model::Model{FT}) where {FT<:AbstractFloat}
    model.variables.waterflux.Runon.RunonUrban = urban_average(
        model.variables.waterflux.Runon.RunonRoofTot,
        model.variables.waterflux.Runon.RunonGroundTot,
        model.parameters.urbangeometry,
    )

    model.variables.waterflux.Runon.RunoffUrban = urban_average(
        model.variables.waterflux.Runon.RunoffRoofTot,
        model.variables.waterflux.Runon.RunoffGroundTot,
        model.parameters.urbangeometry,
    )

    model.variables.radiationflux.SWRabs.TotalUrban = urban_average(
        model.variables.radiationflux.SWRabs.TotalRoof,
        model.variables.radiationflux.SWRabs.TotalCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.radiationflux.SWRin.TotalUrban = urban_average(
        model.variables.radiationflux.SWRin.TotalRoof,
        model.variables.radiationflux.SWRin.TotalCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.radiationflux.SWRout.TotalUrban = urban_average(
        model.variables.radiationflux.SWRout.TotalRoof,
        model.variables.radiationflux.SWRout.TotalCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.radiationflux.LWRabs.TotalUrban = urban_average(
        model.variables.radiationflux.LWRabs.TotalRoof,
        model.variables.radiationflux.LWRabs.TotalCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.radiationflux.LWRin.TotalUrban = urban_average(
        model.variables.radiationflux.LWRin.TotalRoof,
        model.variables.radiationflux.LWRin.TotalCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.radiationflux.LWRout.TotalUrban = urban_average(
        model.variables.radiationflux.LWRout.TotalRoof,
        model.variables.radiationflux.LWRout.TotalCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.heatflux.Hflux.HfluxUrban = urban_average(
        model.variables.heatflux.Hflux.HfluxRoof,
        model.variables.heatflux.Hflux.HfluxCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.heatflux.LEflux.LEfluxUrban = urban_average(
        model.variables.heatflux.LEflux.LEfluxRoof,
        model.variables.heatflux.LEflux.LEfluxCanyon,
        model.parameters.urbangeometry,
    )

    model.variables.heatflux.Gflux.G1Urban = urban_average(
        model.variables.heatflux.Gflux.G1Roof,
        model.variables.heatflux.Gflux.G1Canyon,
        model.parameters.urbangeometry,
    )

    model.variables.heatflux.Gflux.G2Urban = urban_average(
        model.variables.heatflux.Gflux.G2Roof,
        model.variables.heatflux.Gflux.G2Canyon,
        model.parameters.urbangeometry,
    )

    model.variables.waterflux.Leakage.LkUrban = urban_average(
        model.variables.waterflux.Leakage.LkRoof,
        model.variables.waterflux.Leakage.LkGround,
        model.parameters.urbangeometry,
    )

    return nothing
end
