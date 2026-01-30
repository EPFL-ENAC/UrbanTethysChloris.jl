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
        model.variables.waterflux.Owater.ORoofSoilVeg[:] .= 0
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

    TempVec_ittm = deepcopy(model.variables.temperature.tempvec)
    Humidity_ittm = deepcopy(model.variables.humidity.Humidity)
    TempVecB_ittm = deepcopy(model.variables.buildingenergymodel.TempVecB)

    TempVec_ittm2Ext = ExtrapolatedTempVec(model.variables.temperature.tempvec)
    Humidity_ittm2Ext = ExtrapolatedHumidity(model.variables.humidity.Humidity)
    TempVecB_ittm2Ext = ExtrapolatedTempVecB(model.variables.buildingenergymodel.TempVecB)
    Meteo_ittm = Meteotm1(model.forcing.meteorological)

    SWRout_t = Radiation.RadiationFluxes(FT)
    LWRout_t = Radiation.RadiationFluxes(FT)
    T2m = FT(NaN)
    RH_T2m = FT(NaN)

    TempDamp_t = deepcopy(model.variables.temperature.tempdamp)
    RES_t = deepcopy(model.variables.environmentalconditions.resistance)
    Interception_t = deepcopy(model.variables.waterflux.Interception)
    ExWater_t = deepcopy(model.variables.waterflux.ExWater)
    Owater_t = deepcopy(model.variables.waterflux.Owater)
    SoilPotW_t = deepcopy(model.variables.waterflux.SoilPotW)
    CiCO2Leaf_t = deepcopy(model.variables.waterflux.CiCO2Leaf)
    Runon_t = deepcopy(model.variables.waterflux.Runon)
    Qinlat_t = deepcopy(model.variables.waterflux.Qinlat)
    Vwater_t = deepcopy(model.variables.waterflux.Vwater)

    results = create_results_struct(
        FT, NN, model.parameters.soil.roof.ms, model.parameters.soil.ground.ms
    )
    results["ViewFactor"] = ViewFactor
    results["OwaterInitial"] = OwaterInitial

    results_dict, accessors = prepare_results(typeof(output_level), model, NN)

    for i in 1:NN
        @info "Starting iteration $i / $NN"

        if i > 1
            update!(Meteo_ittm, model.forcing.meteorological)
            extrapolate!(TempVec_ittm2Ext, model.variables.temperature.tempvec, i)
            extrapolate!(Humidity_ittm2Ext, model.variables.humidity.Humidity, i)
            extrapolate!(TempVecB_ittm2Ext, model.variables.buildingenergymodel.TempVecB, i)
        end

        if RESPreCalc || fconvPreCalc
            fconv, rsRoofPreCalc, rsGroundPreCalc, rsTreePreCalc = Resistance.precalculate_for_faster_numerical_solution(
                model, i, 1, ViewFactor, BEM_on
            )
        else
            fconv = FT(NaN)
            rsRoofPreCalc = (;)
            rsGroundPreCalc = (;)
            rsTreePreCalc = (;)
        end

        ParHVAC, ParHVACorig = BuildingEnergyModel.ac_heating_turn_on_off(model, BEM_on)

        EnergyUse = (;);
        Fluxes = nothing
        WaterFluxes = nothing
        LEbuildIntc = (;)
        HbuildIntc = (;)
        GbuildIntc = (;)
        WasteHeat = (;)

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
            SWRabsRoofImp, SWRabsRoofVeg, SWRabsTotalRoof, SWRoutRoofImp, SWRoutRoofVeg, SWRoutTotalRoof, SWRinRoofImp, SWRinRoofVeg, SWRinTotalRoof, SWREBRoofImp, SWREBRoofVeg, SWREBTotalRoof, LWRabsRoofVeg, LWRabsRoofImp, LWRabsTotalRoof, LWRoutRoofVeg, LWRoutRoofImp, LWRoutTotalRoof, LWRinRoofImp, LWRinRoofVeg, LWRinTotalRoof, LWREBRoofImp, LWREBRoofVeg, LWREBTotalRoof, HfluxRoofImp, HfluxRoofVeg, HfluxRoof, LEfluxRoofImp, LEfluxRoofVegInt, LEfluxRoofVegPond, LEfluxRoofVegSoil, LTEfluxRoofVeg, LEfluxRoofVeg, LEfluxRoof, G1RoofImp, G2RoofImp, dsRoofImp, G1RoofVeg, G2RoofVeg, dsRoofVeg, G1Roof, G2Roof, dsRoof, raRooftoAtm, rb_LRoof, rap_LRoof, r_soilRoof, rs_sunRoof, rs_shdRoof, EfluxRoofImp, EfluxRoofVegInt, EfluxRoofVegPond, EfluxRoofVegSoil, TEfluxRoofVeg, EfluxRoofVeg, EfluxRoof, QRoofImp, QRoofVegDrip, QRoofVegPond, LkRoofImp, LkRoofVeg, LkRoof, QRoofVegSoil, RunoffRoofTot, RunonRoofTot, IntRoofImp, IntRoofVegPlant, IntRoofVegGround, dInt_dtRoofImp, dInt_dtRoofVegPlant, dInt_dtRoofVegGround, IntRooftot, dInt_dtRooftot, dVRoofSoilVeg_dt, fRoofVeg, VRoofSoilVeg, OwRoofSoilVeg, OSwRoofSoilVeg, ExWaterRoofVeg_H, SoilPotWRoofVeg_H, SoilPotWRoofVeg_L, ExWaterRoofVeg_L, CiCO2LeafRoofVegSun, CiCO2LeafRoofVegShd, WBRoofVegInVeg, WBRoofVegInGround, WBRoofVegSoil, EBRoofImp, EBRoofVeg, Yroof, WBRoofImp, WBRoofVeg, WBRoofTot = eb_wb_roof!(
                model,
                TR,
                TB,
                TempVec_ittm,
                ParCalculation,
                BEM_on,
                RESPreCalc,
                rsRoofPreCalc,
            )

            SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t, albedo_canyon, LWRin_t, LWRout_t, LWRabs_t, LWREB_t, HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, HfluxGround, EfluxGroundImp, EfluxGroundBarePond, EfluxGroundBareSoil, EfluxGroundVegInt, EfluxGroundVegPond, EfluxGroundVegSoil, TEfluxGroundVeg, EfluxTreeInt, TEfluxTree, EfluxGroundBare, EfluxGroundVeg, EfluxGround, EfluxTree, LEfluxGroundImp, LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg, LEfluxTreeInt, LTEfluxTree, LEfluxGroundBare, LEfluxGroundVeg, LEfluxGround, LEfluxTree, CiCO2LeafTreeSun, CiCO2LeafTreeShd, CiCO2LeafGroundVegSun, CiCO2LeafGroundVegShd, raCanyontoAtm, raCanyontoAtmOrig, rap_can, rap_Htree_In, rb_HGround, rb_LGround, r_soilGroundbare, r_soilGroundveg, alp_soilGroundbare, alp_soilGroundveg, rs_sunGround, rs_shdGround, rs_sunTree, rs_shdTree, Fsun_L, Fshd_L, dw_L, RES_w1, RES_w2, rap_W1_In, rap_W2_In, rap_Zp1, HfluxWallSun, HfluxWallShade, EfluxWallSun, EfluxWallShade, LEfluxWallSun, LEfluxWallShade, HfluxCanyon, LEfluxCanyon, EfluxCanyon, G1WallSun, G2WallSun, dsWallSun, G1WallShade, G2WallShade, dsWallShade, G1GroundImp, TDampGroundImp, G1GroundBare, TDampGroundBare, G1GroundVeg, TDampGroundVeg, GTree, TDampTree, G1Ground, G1Canyon, G2Canyon, dsGroundImp, dsGroundBare, dsGroundVeg, dsTree, dsCanyonAir, Ycanyon, QTree, IntTree, dInt_dtTree, QGroundVegDrip, IntGroundVegPlant, dInt_dtGroundVegPlant, QGroundImp, IntGroundImp, dInt_dtGroundImp, fGroundImp, QGroundBarePond, IntGroundBare, dInt_dtGroundBare, fGroundBare, QGroundVegPond, IntGroundVegGround, dInt_dtGroundVegGround, fGroundVeg, VGroundSoilImp, OwGroundSoilImp, OSwGroundSoilImp, LkGroundImp, SoilPotWGroundImp_H, SoilPotWGroundImp_L, ExWaterGroundImp_H, ExWaterGroundImp_L, Rd_gimp, TEgveg_imp, TEtree_imp, Egimp_soil, dVGroundSoilImp_dt, Psi_Soil_gimp, Kf_gimp, VGroundSoilBare, OwGroundSoilBare, OSwGroundSoilBare, LkGroundBare, SoilPotWGroundBare_H, SoilPotWGroundBare_L, ExWaterGroundBare_H, ExWaterGroundBare_L, QGroundBareSoil, TEgveg_bare, TEtree_bare, Egbare_Soil, dVGroundSoilBare_dt, Psi_soil_gbare, Kf_gbare, VGroundSoilVeg, OwGroundSoilVeg, OSwGroundSoilVeg, LkGroundVeg, SoilPotWGroundVeg_H, SoilPotWGroundVeg_L, ExWaterGroundVeg_H, ExWaterGroundVeg_L, QGroundVegSoil, TEgveg_veg, TEtree_veg, Egveg_Soil, dVGroundSoilVeg_dt, Psi_soil_gveg, Kf_gveg, Qin_imp, Qin_bare, Qin_veg, Qin_bare2imp, Qin_bare2veg, Qin_imp2bare, Qin_imp2veg, Qin_veg2imp, Qin_veg2bare, VGroundSoilTot, OwGroundSoilTot, OSwGroundSoilTot, LkGround, Rd, dVGroundSoilTot_dt, SoilPotWGroundTot_L, ExWaterGroundTot_L, TEgveg_tot, SoilPotWGroundTot_H, ExWaterGroundTot_H, TEtree_tot, EB_TEtree, EB_TEgveg, WBIndv, WBTot, RunoffGroundTot, RunonGroundTot, Etot, DeepGLk, StorageTot, EBGroundImp, EBGroundBare, EBGroundVeg, EBTree, EBWallSun, EBWallShade, EBWallSunInt, EBWallShadeInt, EBCanyonT, EBCanyonQ, HumidityCan, HumidityAtm, u_Hcan, u_Zref_und, T2m, q2m, e_T2m, RH_T2m, qcan, e_Tcan, RH_Tcan, DHi, Himp_2m, Hbare_2m, Hveg_2m, Hwsun_2m, Hwshade_2m, Hcan_2m, DEi, Eimp_2m, Ebare_soil_2m, Eveg_int_2m, Eveg_soil_2m, TEveg_2m, Ecan_2m, dS_H_air, dS_LE_air = eb_wb_canyon!(
                model,
                TC,
                TB,
                TempVec_ittm,
                Humidity_ittm,
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

            HbuildIntc, LEbuildIntc, GbuildIntc, SWRabsBc, LWRabsBc, TDampGroundBuild, WasteHeat, EnergyUse, ParACHeat_t, YBuildInt = BuildingEnergyModel.eb_solver_building_output(
                model,
                TC,
                TB,
                TempVecB_ittm,
                TempVec_ittm,
                Humidity_ittm,
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

            # TODO: check why we need Humidity_t, and not directly updating the model?
            Humidity_t.CanyonRelative = HumidityCan.CanyonRelative
            Humidity_t.CanyonSpecific = HumidityCan.CanyonSpecific
            Humidity_t.CanyonVapourPre = HumidityCan.CanyonVapourPre
            Humidity_t.CanyonRelativeSat = HumidityCan.CanyonRelativeSat
            Humidity_t.CanyonSpecificSat = HumidityCan.CanyonSpecificSat
            Humidity_t.CanyonVapourPreSat = HumidityCan.CanyonVapourPreSat
            Humidity_t.AtmRelative = model.forcing.meteorological.AtmRelative
            Humidity_t.AtmSpecific = model.forcing.meteorological.AtmSpecific
            Humidity_t.AtmVapourPre = model.forcing.meteorological.AtmVapourPre
            Humidity_t.AtmRelativeSat = model.forcing.meteorological.AtmRelativeSat
            Humidity_t.AtmSpecificSat = model.forcing.meteorological.AtmSpecificSat
            Humidity_t.AtmVapourPreSat = model.forcing.meteorological.AtmVapourPreSat

            # TODO: modify the structure in-place
            TempDamp_t.TDampGroundImp = TDampGroundImp
            TempDamp_t.TDampGroundBare = TDampGroundBare
            TempDamp_t.TDampGroundVeg = TDampGroundVeg
            TempDamp_t.TDampTree = TDampTree
            TempDamp_t.TDampGroundBuild = TDampGroundBuild

            # RES
            RES_t.raRooftoAtm = raRooftoAtm
            RES_t.raCanyontoAtm = raCanyontoAtm
            RES_t.rap_LRoof = rap_LRoof
            RES_t.rb_LRoof = rb_LRoof
            RES_t.r_soilRoof = r_soilRoof
            RES_t.rs_sunRoof = rs_sunRoof
            RES_t.rs_shdRoof = rs_shdRoof
            RES_t.raCanyontoAtmOrig = raCanyontoAtmOrig
            RES_t.rap_can = rap_can
            RES_t.rap_Htree_In = rap_Htree_In
            RES_t.rb_HGround = rb_HGround
            RES_t.rb_LGround = rb_LGround
            RES_t.r_soilGroundbare = r_soilGroundbare
            RES_t.r_soilGroundveg = r_soilGroundveg
            RES_t.alp_soilGroundbare = alp_soilGroundbare
            RES_t.alp_soilGroundveg = alp_soilGroundveg
            RES_t.rs_sunGround = rs_sunGround
            RES_t.rs_shdGround = rs_shdGround
            RES_t.rs_sunTree = rs_sunTree
            RES_t.rs_shdTree = rs_shdTree
            RES_t.RES_w1 = RES_w1
            RES_t.RES_w2 = RES_w2
            RES_t.rap_W1_In = rap_W1_In
            RES_t.rap_W2_In = rap_W2_In
            RES_t.rap_Zp1 = rap_Zp1

            # Interception
            Interception_t.IntRoofImp = IntRoofImp
            Interception_t.IntRoofVegPlant = IntRoofVegPlant
            Interception_t.IntRoofVegGround = IntRoofVegGround
            Interception_t.IntRooftot = IntRooftot
            Interception_t.IntGroundVegPlant = IntGroundVegPlant
            Interception_t.IntGroundImp = IntGroundImp
            Interception_t.IntGroundBare = IntGroundBare
            Interception_t.IntGroundVegGround = IntGroundVegGround
            Interception_t.IntTree = IntTree

            # ExWater
            ExWater_t.ExWaterRoofVeg_H[:] = ExWaterRoofVeg_H
            ExWater_t.ExWaterRoofVeg_L[:] = ExWaterRoofVeg_L
            ExWater_t.ExWaterGroundImp_H[:] = ExWaterGroundImp_H
            ExWater_t.ExWaterGroundImp_L[:] = ExWaterGroundImp_L
            ExWater_t.ExWaterGroundBare_H[:] = ExWaterGroundBare_H
            ExWater_t.ExWaterGroundBare_L[:] = ExWaterGroundBare_L
            ExWater_t.ExWaterGroundVeg_H[:] = ExWaterGroundVeg_H
            ExWater_t.ExWaterGroundVeg_L[:] = ExWaterGroundVeg_L
            ExWater_t.ExWaterGroundTot_H[:] = ExWaterGroundTot_H
            ExWater_t.ExWaterGroundTot_L[:] = ExWaterGroundTot_L

            # Vwater
            Vwater_t.VRoofSoilVeg[:] = VRoofSoilVeg
            Vwater_t.VGroundSoilImp[:] = VGroundSoilImp
            Vwater_t.VGroundSoilBare[:] = VGroundSoilBare
            Vwater_t.VGroundSoilVeg[:] = VGroundSoilVeg
            Vwater_t.VGroundSoilTot[:] = VGroundSoilTot

            # Owater
            Owater_t.OwRoofSoilVeg[:] = OwRoofSoilVeg
            Owater_t.OwGroundSoilImp[:] = OwGroundSoilImp
            Owater_t.OwGroundSoilBare[:] = OwGroundSoilBare
            Owater_t.OwGroundSoilVeg[:] = OwGroundSoilVeg
            Owater_t.OwGroundSoilTot[:] = OwGroundSoilTot

            # SoilPotW
            SoilPotW_t.SoilPotWRoofVeg_H = SoilPotWRoofVeg_H
            SoilPotW_t.SoilPotWRoofVeg_L = SoilPotWRoofVeg_L
            SoilPotW_t.SoilPotWGroundImp_H = SoilPotWGroundImp_H
            SoilPotW_t.SoilPotWGroundImp_L = SoilPotWGroundImp_L
            SoilPotW_t.SoilPotWGroundBare_H = SoilPotWGroundBare_H
            SoilPotW_t.SoilPotWGroundBare_L = SoilPotWGroundBare_L
            SoilPotW_t.SoilPotWGroundVeg_H = SoilPotWGroundVeg_H
            SoilPotW_t.SoilPotWGroundVeg_L = SoilPotWGroundVeg_L
            SoilPotW_t.SoilPotWGroundTot_L = SoilPotWGroundTot_L
            SoilPotW_t.SoilPotWGroundTot_H = SoilPotWGroundTot_H

            # CiCO2Leaf
            CiCO2Leaf_t.CiCO2LeafRoofVegSun = CiCO2LeafRoofVegSun
            CiCO2Leaf_t.CiCO2LeafRoofVegShd = CiCO2LeafRoofVegShd
            CiCO2Leaf_t.CiCO2LeafGroundVegSun = CiCO2LeafGroundVegSun
            CiCO2Leaf_t.CiCO2LeafGroundVegShd = CiCO2LeafGroundVegShd
            CiCO2Leaf_t.CiCO2LeafTreeSun = CiCO2LeafTreeSun
            CiCO2Leaf_t.CiCO2LeafTreeShd = CiCO2LeafTreeShd

            # Runon
            Runon_t.RunonRoofTot = RunonRoofTot
            Runon_t.RunoffRoofTot = RunoffRoofTot
            Runon_t.RunonGroundTot = RunonGroundTot
            Runon_t.RunoffGroundTot = RunoffGroundTot
            Runon_t.RunonUrban = urban_average(
                RunonRoofTot, RunonGroundTot, model.parameters.urbangeometry
            )
            Runon_t.RunoffUrban = urban_average(
                RunoffRoofTot, RunoffGroundTot, model.parameters.urbangeometry
            )

            # Qinlat
            Qinlat_t.Qin_bare2imp = Qin_bare2imp
            Qinlat_t.Qin_veg2imp = Qin_veg2imp
            Qinlat_t.Qin_veg2bare = Qin_veg2bare
            Qinlat_t.Qin_imp2bare = Qin_imp2bare
            Qinlat_t.Qin_bare2veg = Qin_bare2veg
            Qinlat_t.Qin_imp2veg = Qin_imp2veg
            Qinlat_t.Qin_imp = Qin_imp
            Qinlat_t.Qin_bare = Qin_bare
            Qinlat_t.Qin_veg = Qin_veg

            Fluxes = (;
                SWRin_t,
                SWRout_t,
                SWRabs_t,
                LWRin_t,
                LWRout_t,
                LWRabs_t,
                SWRinTotalRoof,
                SWRabsTotalRoof,
                SWRoutTotalRoof,
                LWRinTotalRoof,
                LWRabsTotalRoof,
                LWRoutTotalRoof,
                SWRabsTotalUrban=urban_average(
                    SWRabsTotalRoof, SWRabs_t.TotalCanyon, model.parameters.urbangeometry
                ),
                SWRinTotalUrban=urban_average(
                    SWRinTotalRoof, SWRin_t.TotalCanyon, model.parameters.urbangeometry
                ),
                SWRoutTotalUrban=urban_average(
                    SWRoutTotalRoof, SWRout_t.TotalCanyon, model.parameters.urbangeometry
                ),
                LWRabsTotalUrban=urban_average(
                    LWRabsTotalRoof, LWRabs_t.TotalCanyon, model.parameters.urbangeometry
                ),
                LWRinTotalUrban=urban_average(
                    LWRinTotalRoof, LWRin_t.TotalCanyon, model.parameters.urbangeometry
                ),
                LWRoutTotalUrban=urban_average(
                    LWRoutTotalRoof, LWRout_t.TotalCanyon, model.parameters.urbangeometry
                ),
                LEfluxRoof,
                HfluxRoof,
                G1Roof,
                G2Roof,
                dsRoof,
                LEfluxCanyon,
                HfluxCanyon,
                G1Canyon,
                G2Canyon,
                G1Ground,
                G1WallSun,
                G1WallShade,
                dsWallSun,
                dsWallShade,
                dS_H_air,
                dS_LE_air,
                HfluxUrban=urban_average(
                    HfluxRoof, HfluxCanyon, model.parameters.urbangeometry
                ),
                LEfluxUrban=urban_average(
                    LEfluxRoof, LEfluxCanyon, model.parameters.urbangeometry
                ),
                G1Urban=urban_average(G1Roof, G1Canyon, model.parameters.urbangeometry),
                G2Urban=urban_average(G2Roof, G2Canyon, model.parameters.urbangeometry),
            )

            # Store water balance components
            Runon_nt = (;
                RunoffRoofTot,
                RunoffGroundTot,
                RunoffUrban=Runon_t.RunoffUrban,
                RunonRoofTot,
                RunonGroundTot,
                RunonUrban=Runon_t.RunonUrban,
            )

            Leakage_nt = (;
                LkRoof,
                LkGround,
                LkUrban=urban_average(LkRoof, LkGround, model.parameters.urbangeometry),
            )

            LEflux_nt = (;
                LEfluxRoofImp,
                LEfluxRoofVegInt,
                LEfluxRoofVegPond,
                LEfluxRoofVegSoil,
                LTEfluxRoofVeg,
                LEfluxGroundImp,
                LEfluxGroundBarePond,
                LEfluxGroundBareSoil,
                LEfluxGroundVegInt,
                LEfluxGroundVegPond,
                LEfluxGroundVegSoil,
                LTEfluxGroundVeg,
                LEfluxTreeInt,
                LTEfluxTree,
            )

            dVwater_dt_nt = (; dVRoofSoilVeg_dt, dVGroundSoilTot_dt)

            dInt_dt_nt = (;
                dInt_dtRoofVegPlant,
                dInt_dtGroundVegPlant,
                dInt_dtTree,
                dInt_dtRoofVegGround,
                dInt_dtRoofImp,
                dInt_dtGroundVegGround,
                dInt_dtGroundBare,
                dInt_dtGroundImp,
            )

            Int_nt = (;
                IntRooftot,
                IntGroundImp,
                IntGroundBare,
                IntGroundVegPlant,
                IntGroundVegGround,
                IntTree,
            )

            LEbuildInt_nt = (;
                LEpeople=LEbuildIntc.LEpeople,
                LEequip=LEbuildIntc.LEequip,
                LEvent=LEbuildIntc.LEvent,
            )

            WaterFluxes = (
                Runon=Runon_nt,
                Leakage=Leakage_nt,
                LEflux=LEflux_nt,
                dVwater_dt=dVwater_dt_nt,
                dInt_dt=dInt_dt_nt,
                Int=Int_nt,
                LEbuildInt=LEbuildInt_nt,
            )
        end

        Tmrt = MeanRadiantTemperature.mean_radiant_temperature!(
            model, SWRout_t, LWRout_t, ViewFactorPoint
        )
        # TODO: check whether we should be using the hour as a float (e.g. 10.5 for 10:30) or
        # the hour as an integer (10 for 10:30).

        # TODO: Modify model in place
        u_ZPerson = Resistance.wind_profile_point_output(model)

        # TODO: modify in place
        model.variables.temperature.thermalcomfort.UTCI = OutdoorThermalComfort.utci_approx(
            T2m - FT(273.15), RH_T2m * 100, Tmrt, u_ZPerson
        )

        # Assign outputs
        model.variables.energybalance.Solver.YfunctionOutput = vcat(
            Yroof, Ycanyon, YBuildInt
        )
        # TODO: add missing TotalRoof field to RadiationFluxes struct
        # TODO: add urban average for Hflux, LEflux, Gflux
        # TODO: implement energy balance check script as function

        # tempvec - already done
        model.variables.temperature.tempvec.T2m = T2m

        # tempdamp
        update!(model.variables.temperature.tempdamp, TempDamp_t)

        # Humidity - careful, update! already exists
        update!(model.variables.humidity.Humidity, Humidity_t)

        # TempVecB - already done higher

        # Resistance
        update!(model.variables.environmentalconditions.resistance, RES_t)

        # Water fluxes
        update!(model.variables.waterflux.Interception, Interception_t)
        update!(model.variables.waterflux.ExWater, ExWater_t)
        update!(model.variables.waterflux.Vwater, Vwater_t)
        update!(
            model.variables.waterflux.Owater,
            Owater_t,
            model.parameters.soil.roof,
            model.parameters.soil.ground,
            O33,
        )
        update!(model.variables.waterflux.SoilPotW, SoilPotW_t)
        update!(model.variables.waterflux.CiCO2Leaf, CiCO2Leaf_t)
        update!(model.variables.waterflux.Runon, Runon_t)
        update!(model.variables.waterflux.Qinlat, Qinlat_t)

        update!(TempVecB_ittm, model.variables.buildingenergymodel.TempVecB)
        update!(TempVec_ittm, model.variables.temperature.tempvec)
        update!(Humidity_ittm, model.variables.humidity.Humidity)

        store_results!(results, model, i)
        store_fluxes!(results, i, Fluxes)
        store_GbuildInt!(results, i, GbuildIntc)
        store_HbuildInt!(results, i, HbuildIntc)
        store_LEbuildInt!(results, i, LEbuildIntc)
        store_BEMWasteHeat!(results, i, WasteHeat)

        store_water_fluxes!(results, i, WaterFluxes)

        assign_results!(results_dict, accessors, model, i)

        # Update forcing parameters for the next step
        model.forcing = forcing[i + 1]
    end

    return results, results_dict
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

function create_results_struct(
    ::Type{FT}, NN::Signed, MR::Signed, MG::Signed
) where {FT<:AbstractFloat}
    results = Dict{String,Any}(
        "TRoofImp" => zeros(FT, NN),
        "TRoofVeg" => zeros(FT, NN),
        "TRoofIntImp" => zeros(FT, NN),
        "TRoofIntVeg" => zeros(FT, NN),
        "TGroundImp" => zeros(FT, NN),
        "TGroundBare" => zeros(FT, NN),
        "TGroundVeg" => zeros(FT, NN),
        "TTree" => zeros(FT, NN),
        "TWallSun" => zeros(FT, NN),
        "TWallShade" => zeros(FT, NN),
        "TWallIntSun" => zeros(FT, NN),
        "TWallIntShade" => zeros(FT, NN),
        "TCanyon" => zeros(FT, NN),
        "Tatm" => zeros(FT, NN),
        "T2m" => zeros(FT, NN),
        "RH_T2m" => zeros(FT, NN),
        "Tmrt" => zeros(FT, NN),
        "RHbin" => zeros(FT, NN),
        "qbin" => zeros(FT, NN),
        "UTCI" => zeros(FT, NN),
        "Tbin" => zeros(FT, NN),
        "ViewFactor" => nothing,
        "SWRinGroundVeg" => zeros(FT, NN),
        "SWRinGroundBare" => zeros(FT, NN),
        "SWRinGroundImp" => zeros(FT, NN),
        "SWRinWallShade" => zeros(FT, NN),
        "SWRinWallSun" => zeros(FT, NN),
        "SWRinTree" => zeros(FT, NN),
        "SWRinTotalUrban" => zeros(FT, NN),
        "SWRabsGroundVeg" => zeros(FT, NN),
        "SWRabsGroundBare" => zeros(FT, NN),
        "SWRabsGroundImp" => zeros(FT, NN),
        "SWRabsWallShade" => zeros(FT, NN),
        "SWRabsWallSun" => zeros(FT, NN),
        "SWRabsTree" => zeros(FT, NN),
        "SWRabsTotalUrban" => zeros(FT, NN),
        "SWRoutGroundVeg" => zeros(FT, NN),
        "SWRoutGroundBare" => zeros(FT, NN),
        "SWRoutGroundImp" => zeros(FT, NN),
        "SWRoutWallShade" => zeros(FT, NN),
        "SWRoutWallSun" => zeros(FT, NN),
        "SWRoutTree" => zeros(FT, NN),
        "SWRoutTotalUrban" => zeros(FT, NN),
        "LWRinGroundVeg" => zeros(FT, NN),
        "LWRinGroundBare" => zeros(FT, NN),
        "LWRinGroundImp" => zeros(FT, NN),
        "LWRinWallShade" => zeros(FT, NN),
        "LWRinWallSun" => zeros(FT, NN),
        "LWRinTree" => zeros(FT, NN),
        "LWRinTotalUrban" => zeros(FT, NN),
        "LWRabsGroundVeg" => zeros(FT, NN),
        "LWRabsGroundBare" => zeros(FT, NN),
        "LWRabsGroundImp" => zeros(FT, NN),
        "LWRabsWallShade" => zeros(FT, NN),
        "LWRabsWallSun" => zeros(FT, NN),
        "LWRabsTree" => zeros(FT, NN),
        "LWRabsTotalUrban" => zeros(FT, NN),
        "LWRoutGroundVeg" => zeros(FT, NN),
        "LWRoutGroundBare" => zeros(FT, NN),
        "LWRoutGroundImp" => zeros(FT, NN),
        "LWRoutWallShade" => zeros(FT, NN),
        "LWRoutWallSun" => zeros(FT, NN),
        "LWRoutTree" => zeros(FT, NN),
        "LWRoutTotalUrban" => zeros(FT, NN),
        "SWRinTotalRoof" => zeros(FT, NN),
        "SWRabsTotalRoof" => zeros(FT, NN),
        "SWRoutTotalRoof" => zeros(FT, NN),
        "LWRinTotalRoof" => zeros(FT, NN),
        "LWRabsTotalRoof" => zeros(FT, NN),
        "LWRoutTotalRoof" => zeros(FT, NN),
        "LEfluxRoof" => zeros(FT, NN),
        "HfluxRoof" => zeros(FT, NN),
        "G1Roof" => zeros(FT, NN),
        "G2Roof" => zeros(FT, NN),
        "dsRoof" => zeros(FT, NN),
        "LEfluxCanyon" => zeros(FT, NN),
        "HfluxCanyon" => zeros(FT, NN),
        "HfluxUrban" => zeros(FT, NN),
        "LEfluxUrban" => zeros(FT, NN),
        "G1Urban" => zeros(FT, NN),
        "G2Urban" => zeros(FT, NN),
        "G1Canyon" => zeros(FT, NN),
        "G2Canyon" => zeros(FT, NN),
        "G1Ground" => zeros(FT, NN),
        "G1WallSun" => zeros(FT, NN),
        "G1WallShade" => zeros(FT, NN),
        "dsWallSun" => zeros(FT, NN),
        "dsWallShade" => zeros(FT, NN),
        "dS_H_air" => zeros(FT, NN),
        "dS_LE_air" => zeros(FT, NN),
        "SWRabsWallSunTransmitted" => zeros(FT, NN),
        "SWRabsWallShadeTransmitted" => zeros(FT, NN),
        "SWRabsTotalCanyon" => zeros(FT, NN),
        "LWRabsTotalCanyon" => zeros(FT, NN),
        "HbuildIntdSH_air" => zeros(FT, NN),
        "LEbuildIntdSLE_air" => zeros(FT, NN),
        "GbuildIntGfloor" => zeros(FT, NN),
        "GbuildIntdSinternalMass" => zeros(FT, NN),
        "WasteHeatTotAnthInput_URB" => zeros(FT, NN),
        "WasteHeatWaterFromAC_Can" => zeros(FT, NN),
        "WasteHeatSensibleFromVent_Can" => zeros(FT, NN),
        "WasteHeatSensibleFromAC_Can" => zeros(FT, NN),
        "WasteHeatSensibleFromHeat_Can" => zeros(FT, NN),
        "WasteHeatLatentFromVent_Can" => zeros(FT, NN),
        "WasteHeatLatentFromAC_Can" => zeros(FT, NN),
        "WasteHeatLatentFromHeat_Can" => zeros(FT, NN),
        # Water balance components - Runoff and leakage
        "RunoffRoofTot" => zeros(FT, NN),
        "RunoffGroundTot" => zeros(FT, NN),
        "RunoffUrban" => zeros(FT, NN),
        "LkRoof" => zeros(FT, NN),
        "LkGround" => zeros(FT, NN),
        "LkUrban" => zeros(FT, NN),
        # Water balance components - LE fluxes by source
        "LEfluxRoofImp" => zeros(FT, NN),
        "LEfluxRoofVegInt" => zeros(FT, NN),
        "LEfluxRoofVegPond" => zeros(FT, NN),
        "LEfluxRoofVegSoil" => zeros(FT, NN),
        "LTEfluxRoofVeg" => zeros(FT, NN),
        "LEfluxGroundImp" => zeros(FT, NN),
        "LEfluxGroundBarePond" => zeros(FT, NN),
        "LEfluxGroundBareSoil" => zeros(FT, NN),
        "LEfluxGroundVegInt" => zeros(FT, NN),
        "LEfluxGroundVegPond" => zeros(FT, NN),
        "LEfluxGroundVegSoil" => zeros(FT, NN),
        "LTEfluxGroundVeg" => zeros(FT, NN),
        "LEfluxTreeInt" => zeros(FT, NN),
        "LTEfluxTree" => zeros(FT, NN),
        # Water balance components - Building LE
        "LEbuildIntLEpeople" => zeros(FT, NN),
        "LEbuildIntLEequip" => zeros(FT, NN),
        "LEbuildIntLEvent" => zeros(FT, NN),
        # Water balance components - Soil moisture change
        "dVRoofSoilVeg_dt" => zeros(FT, NN),
        "dVGroundSoilTot_dt" => zeros(FT, NN),
        # Water balance components - Interception change
        "dInt_dtRoofVegPlant" => zeros(FT, NN),
        "dInt_dtGroundVegPlant" => zeros(FT, NN),
        "dInt_dtTree" => zeros(FT, NN),
        "dInt_dtRoofVegGround" => zeros(FT, NN),
        "dInt_dtRoofImp" => zeros(FT, NN),
        "dInt_dtGroundVegGround" => zeros(FT, NN),
        "dInt_dtGroundBare" => zeros(FT, NN),
        "dInt_dtGroundImp" => zeros(FT, NN),
        # Water balance components - Runon
        "RunonRoofTot" => zeros(FT, NN),
        "RunonGroundTot" => zeros(FT, NN),
        "RunonUrban" => zeros(FT, NN),
        # Water balance components - Interception storage
        "IntRooftot" => zeros(FT, NN),
        "IntGroundImp" => zeros(FT, NN),
        "IntGroundBare" => zeros(FT, NN),
        "IntGroundVegPlant" => zeros(FT, NN),
        "IntGroundVegGround" => zeros(FT, NN),
        "IntTree" => zeros(FT, NN),
        "OwRoofSoilVeg" => zeros(FT, MR, NN),
        "OwGroundSoilImp" => zeros(FT, MG, NN),
        "OwGroundSoilBare" => zeros(FT, MG, NN),
        "OwGroundSoilVeg" => zeros(FT, MG, NN),
    )

    return results
end

function store_results!(
    results::Dict{String,Any}, model::Model{FT}, i::Signed
) where {FT<:AbstractFloat}
    store_results!(results, model.variables.temperature.tempvec, i)
    store_results!(results, model.variables.buildingenergymodel.HumidityBuilding, i)
    store_results!(results, model.variables.temperature.mrt, i)
    store_results!(results, model.variables.temperature.thermalcomfort, i)
    store_results!(results, model.variables.humidity.Results2m, i)
    store_results!(results, model.variables.buildingenergymodel.TempVecB, i)
    store_results!(results, model.variables.waterflux.Owater, i)
end

function store_results!(
    results::Dict{String,Any}, Owater::ModelComponents.ModelVariables.Owater{FT}, i::Signed
) where {FT<:AbstractFloat}
    results["OwRoofSoilVeg"][:, i] = Owater.OwRoofSoilVeg
    results["OwGroundSoilImp"][:, i] = Owater.OwGroundSoilImp
    results["OwGroundSoilBare"][:, i] = Owater.OwGroundSoilBare
    results["OwGroundSoilVeg"][:, i] = Owater.OwGroundSoilVeg

    return nothing
end

function store_results!(
    results::Dict{String,Any},
    TempVec::ModelComponents.ModelVariables.TempVec{FT},
    i::Signed,
) where {FT<:AbstractFloat}
    results["TRoofImp"][i] = TempVec.TRoofImp
    results["TRoofVeg"][i] = TempVec.TRoofVeg
    results["TRoofIntImp"][i] = TempVec.TRoofIntImp
    results["TRoofIntVeg"][i] = TempVec.TRoofIntVeg
    results["TGroundImp"][i] = TempVec.TGroundImp
    results["TGroundBare"][i] = TempVec.TGroundBare
    results["TGroundVeg"][i] = TempVec.TGroundVeg
    results["TTree"][i] = TempVec.TTree
    results["TWallSun"][i] = TempVec.TWallSun
    results["TWallShade"][i] = TempVec.TWallShade
    results["TWallIntSun"][i] = TempVec.TWallIntSun
    results["TWallIntShade"][i] = TempVec.TWallIntShade
    results["TCanyon"][i] = TempVec.TCanyon
    results["Tatm"][i] = TempVec.Tatm
    results["T2m"][i] = TempVec.T2m

    return nothing
end

function store_results!(
    results::Dict{String,Any},
    HumidityBuilding::ModelComponents.ModelVariables.HumidityBuilding{FT},
    i::Signed,
) where {FT<:AbstractFloat}
    results["RHbin"][i] = HumidityBuilding.RHbin
    results["qbin"][i] = HumidityBuilding.qbin

    return nothing
end

function store_results!(
    results::Dict{String,Any},
    ThermalComfort::ModelComponents.ModelVariables.ThermalComfort{FT},
    i::Signed,
) where {FT<:AbstractFloat}
    results["UTCI"][i] = ThermalComfort.UTCI

    return nothing
end

function store_results!(
    results::Dict{String,Any}, MRT::ModelComponents.ModelVariables.MRT{FT}, i::Signed
) where {FT<:AbstractFloat}
    results["Tmrt"][i] = MRT.Tmrt

    return nothing
end

function store_results!(
    results::Dict{String,Any},
    TempVecB::ModelComponents.ModelVariables.TempVecB{FT},
    i::Signed,
) where {FT<:AbstractFloat}
    results["Tbin"][i] = TempVecB.Tbin

    return nothing
end

function store_results!(
    results::Dict{String,Any},
    Results2m::ModelComponents.ModelVariables.Results2m{FT},
    i::Signed,
) where {FT<:AbstractFloat}
    results["T2m"][i] = Results2m.T2m
    results["RH_T2m"][i] = Results2m.RH_T2m

    return nothing
end

function store_results!(
    results::Dict{String,Any},
    RadiationFluxes::Radiation.AbstractRadiationFluxes{FT},
    i::Signed,
    prefix::AbstractString,
) where {FT<:AbstractFloat}
    for field in propertynames(RadiationFluxes)
        key_str = prefix * String(field)
        if haskey(results, key_str)
            val = getfield(RadiationFluxes, field)
            results[key_str][i] = val
        end
    end
    return nothing
end

function store_GbuildInt!(results::Dict{String,Any}, i::Signed, GbuildInt::NamedTuple)
    results["GbuildIntGfloor"][i] = GbuildInt.Gfloor
    results["GbuildIntdSinternalMass"][i] = GbuildInt.dSinternalMass
    return nothing
end

function store_HbuildInt!(results::Dict{String,Any}, i::Signed, HbuildInt::NamedTuple)
    results["HbuildIntdSH_air"][i] = HbuildInt.dSH_air
    return nothing
end

function store_LEbuildInt!(results::Dict{String,Any}, i::Signed, LEbuildInt::NamedTuple)
    results["LEbuildIntdSLE_air"][i] = LEbuildInt.dSLE_air
    return nothing
end

# TODO: switch from NamedTuple to BEMWasteHeat composite type!
function store_BEMWasteHeat!(results::Dict{String,Any}, i::Signed, WasteHeat::NamedTuple)
    results["WasteHeatTotAnthInput_URB"][i] = WasteHeat.TotAnthInput_URB
    results["WasteHeatWaterFromAC_Can"][i] = WasteHeat.WaterFromAC_Can
    results["WasteHeatSensibleFromVent_Can"][i] = WasteHeat.SensibleFromVent_Can
    results["WasteHeatSensibleFromAC_Can"][i] = WasteHeat.SensibleFromAC_Can
    results["WasteHeatSensibleFromHeat_Can"][i] = WasteHeat.SensibleFromHeat_Can
    results["WasteHeatLatentFromVent_Can"][i] = WasteHeat.LatentFromVent_Can
    results["WasteHeatLatentFromAC_Can"][i] = WasteHeat.LatentFromAC_Can
    results["WasteHeatLatentFromHeat_Can"][i] = WasteHeat.LatentFromHeat_Can
    return nothing
end

function store_fluxes!(results::Dict{String,Any}, i::Signed, Fluxes::NamedTuple)
    radiation_keys = [:SWRin_t, :SWRout_t, :SWRabs_t, :LWRin_t, :LWRout_t, :LWRabs_t]
    keys_to_store = setdiff(keys(Fluxes), radiation_keys)
    for key in keys_to_store
        val = getproperty(Fluxes, key)
        key_str = String(key)
        if haskey(results, key_str)
            # Assuming the vector is already allocated
            results[key_str][i] = val
        end
    end

    for key in radiation_keys
        RadiationFluxes = getproperty(Fluxes, key)
        prefix = split(String(key), "_")[1]
        store_results!(results, RadiationFluxes, i, prefix)
    end
    return nothing
end

function store_water_fluxes!(results::Dict{String,Any}, i::Signed, WaterFluxes::NamedTuple)
    store_Runon!(results, i, WaterFluxes.Runon)
    store_Leakage!(results, i, WaterFluxes.Leakage)
    store_LEflux!(results, i, WaterFluxes.LEflux)
    store_dVwater_dt!(results, i, WaterFluxes.dVwater_dt)
    store_dInt_dt!(results, i, WaterFluxes.dInt_dt)
    store_Int!(results, i, WaterFluxes.Int)
    store_LEbuildInt_water!(results, i, WaterFluxes.LEbuildInt)

    return nothing
end

function store_LEbuildInt_water!(
    results::Dict{String,Any}, i::Signed, LEbuildInt::NamedTuple
)
    results["LEbuildIntLEpeople"][i] = LEbuildInt.LEpeople
    results["LEbuildIntLEequip"][i] = LEbuildInt.LEequip
    results["LEbuildIntLEvent"][i] = LEbuildInt.LEvent

    return nothing
end

function store_Int!(results::Dict{String,Any}, i::Signed, Int::NamedTuple)
    results["IntRooftot"][i] = Int.IntRooftot
    results["IntGroundImp"][i] = Int.IntGroundImp
    results["IntGroundBare"][i] = Int.IntGroundBare
    results["IntGroundVegPlant"][i] = Int.IntGroundVegPlant
    results["IntGroundVegGround"][i] = Int.IntGroundVegGround
    results["IntTree"][i] = Int.IntTree

    return nothing
end

function store_dInt_dt!(results::Dict{String,Any}, i::Signed, dInt_dt::NamedTuple)
    # Interception change
    dint_keys = [
        :dInt_dtRoofVegPlant,
        :dInt_dtGroundVegPlant,
        :dInt_dtTree,
        :dInt_dtRoofVegGround,
        :dInt_dtRoofImp,
        :dInt_dtGroundVegGround,
        :dInt_dtGroundBare,
        :dInt_dtGroundImp,
    ]

    for key in dint_keys
        results[String(key)][i] = getfield(dInt_dt, key)
    end

    return nothing
end

function store_dVwater_dt!(results::Dict{String,Any}, i::Signed, dVwater_dt::NamedTuple)
    results["dVRoofSoilVeg_dt"][i] = dVwater_dt.dVRoofSoilVeg_dt
    results["dVGroundSoilTot_dt"][i] = dVwater_dt.dVGroundSoilTot_dt

    return nothing
end

function store_LEflux!(results::Dict{String,Any}, i::Signed, LEflux::NamedTuple)
    # LE fluxes by source
    le_flux_keys = [
        :LEfluxRoofImp,
        :LEfluxRoofVegInt,
        :LEfluxRoofVegPond,
        :LEfluxRoofVegSoil,
        :LTEfluxRoofVeg,
        :LEfluxGroundImp,
        :LEfluxGroundBarePond,
        :LEfluxGroundBareSoil,
        :LEfluxGroundVegInt,
        :LEfluxGroundVegPond,
        :LEfluxGroundVegSoil,
        :LTEfluxGroundVeg,
        :LEfluxTreeInt,
        :LTEfluxTree,
    ]
    for key in le_flux_keys
        results[String(key)][i] = getfield(LEflux, key)
    end

    return nothing
end

function store_Leakage!(results::Dict{String,Any}, i::Signed, Leakage::NamedTuple)
    results["LkRoof"][i] = Leakage.LkRoof
    results["LkGround"][i] = Leakage.LkGround
    results["LkUrban"][i] = Leakage.LkUrban

    return nothing
end

function store_Runon!(results::Dict{String,Any}, i::Signed, Runon::NamedTuple)
    results["RunoffRoofTot"][i] = Runon.RunoffRoofTot
    results["RunoffGroundTot"][i] = Runon.RunoffGroundTot
    results["RunoffUrban"][i] = Runon.RunoffUrban
    results["RunonRoofTot"][i] = Runon.RunonRoofTot
    results["RunonGroundTot"][i] = Runon.RunonGroundTot
    results["RunonUrban"][i] = Runon.RunonUrban

    return nothing
end
