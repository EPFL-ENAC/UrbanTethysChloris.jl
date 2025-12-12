function run_simulation(
    model::Model{FT};
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

    Ttot, fval, exitflag = nothing, nothing, nothing

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
    Owater_t = deepcopy(model.variables.waterflux.Owater)

    for i in 1:NN
        @info "Starting iteration $i / $NN"

        if i > 1
            update!(Meteo_ittm, model.forcing.meteorological)
        end

        if i > 2
            extrapolate!(TempVec_ittm2Ext, model.variables.temperature.tempvec)
            extrapolate!(Humidity_ittm2Ext, model.variables.humidity.Humidity)
            extrapolate!(TempVecB_ittm2Ext, model.variables.buildingenergymodel.TempVecB)
        end

        if RESPreCalc || fconvPreCalc
            fconv, rsRoofPreCalc, rsGroundPreCalc, rsTreePreCalc = Resistance.precalculate_for_faster_numerical_solution(
                i,
                1,
                model.variables.temperature.tempvec,
                model.variables.humidity.Humidity,
                model.parameters.vegetation.ground,
                model.variables.waterflux.SoilPotW,
                model.variables.waterflux.CiCO2Leaf,
                model.forcing.meteorological,
                model.parameters.urbangeometry,
                model.parameters.surfacefractions.ground,
                model.parameters.optical.ground,
                model.parameters.optical.wall,
                model.parameters.optical.tree,
                model.parameters.vegetation.tree,
                model.forcing.sunposition,
                ViewFactor,
                model.parameters.building_energy.windows,
                BEM_on,
                model.parameters.vegetation.roof,
                model.parameters.optical.roof,
                model.parameters.surfacefractions.roof,
                model.variables.environmentalconditions.resistance,
            )
        else
            fconv = FT(NaN)
            rsRoofPreCalc = (;)
            rsGroundPreCalc = (;)
            rsTreePreCalc = (;)
        end

        ParHVAC, ParHVACorig = BuildingEnergyModel.ac_heating_turn_on_off(
            model.parameters.building_energy.hvac,
            model.variables.buildingenergymodel.TempVecB,
            model.variables.temperature.tempvec,
            model.variables.humidity.Humidity,
            model.forcing.meteorological,
            model.parameters.urbangeometry,
            BEM_on,
        )

        EnergyUse = (;);

        for HVACittm in 1:2
            if BEM_on && HVACittm == 2
                ParHVAC = update_hvac_parameters(
                    ParHVACorig,
                    ParHVAC,
                    EnergyUse,
                    model.variables.buildingenergymodel.TempVecB.Tbin,
                    model.variables.buildingenergymodel.TempVecB.qbin,
                )
            end
            Ttot, fval, exitflag = f_solver_tot(
                model.variables.temperature.tempvec,
                model.variables.buildingenergymodel.TempVecB,
                model.variables.humidity.Humidity,
                model.forcing.meteorological,
                model.variables.waterflux.Interception,
                model.variables.waterflux.ExWater,
                model.variables.waterflux.Vwater,
                model.variables.waterflux.Owater,
                model.variables.waterflux.SoilPotW,
                model.variables.waterflux.CiCO2Leaf,
                model.variables.temperature.tempdamp,
                ViewFactor,
                model.parameters.urbangeometry,
                model.parameters.surfacefractions.ground,
                model.parameters.surfacefractions.roof,
                WallLayers,
                model.parameters.soil.ground,
                ParInterceptionTree,
                model.parameters.optical.ground,
                model.parameters.optical.wall,
                model.parameters.optical.tree,
                model.parameters.thermal.ground,
                model.parameters.thermal.wall,
                model.parameters.vegetation.ground,
                model.parameters.vegetation.tree,
                model.parameters.soil.roof,
                model.parameters.optical.roof,
                model.parameters.thermal.roof,
                model.parameters.vegetation.roof,
                model.forcing.sunposition,
                model.forcing.meteorological,
                model.forcing.anthropogenic,
                ParCalculation,
                model.parameters.building_energy.indoor_optical,
                ParHVAC,
                model.parameters.building_energy.thermal,
                model.parameters.building_energy.windows,
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
                model.forcing.hvacschedule,
            )

            update!(model.variables.temperature.tempvec, Ttot)
            update!(model.variables.humidity.Humidity, Ttot)
            update!(model.variables.buildingenergymodel.TempVecB, Ttot)

            TR = roof_temperature(model.variables.temperature.tempvec)
            TC = canyon_temperature(
                model.variables.temperature.tempvec, model.variables.humidity.Humidity
            )
            TB = building_temperature(model.variables.buildingenergymodel.TempVecB)

            SWRabsRoofImp, SWRabsRoofVeg, SWRabsTotalRoof, SWRoutRoofImp, SWRoutRoofVeg, SWRoutTotalRoof, SWRinRoofImp, SWRinRoofVeg, SWRinTotalRoof, SWREBRoofImp, SWREBRoofVeg, SWREBTotalRoof, LWRabsRoofVeg, LWRabsRoofImp, LWRabsTotalRoof, LWRoutRoofVeg, LWRoutRoofImp, LWRoutTotalRoof, LWRinRoofImp, LWRinRoofVeg, LWRinTotalRoof, LWREBRoofImp, LWREBRoofVeg, LWREBTotalRoof, HfluxRoofImp, HfluxRoofVeg, HfluxRoof, LEfluxRoofImp, LEfluxRoofVegInt, LEfluxRoofVegPond, LEfluxRoofVegSoil, LTEfluxRoofVeg, LEfluxRoofVeg, LEfluxRoof, G1RoofImp, G2RoofImp, dsRoofImp, G1RoofVeg, G2RoofVeg, dsRoofVeg, G1Roof, G2Roof, dsRoof, raRooftoAtm, rb_LRoof, rap_LRoof, r_soilRoof, rs_sunRoof, rs_shdRoof, EfluxRoofImp, EfluxRoofVegInt, EfluxRoofVegPond, EfluxRoofVegSoil, TEfluxRoofVeg, EfluxRoofVeg, EfluxRoof, QRoofImp, QRoofVegDrip, QRoofVegPond, LkRoofImp, LkRoofVeg, LkRoof, QRoofVegSoil, RunoffRoofTot, RunonRoofTot, IntRoofImp, IntRoofVegPlant, IntRoofVegGround, dInt_dtRoofImp, dInt_dtRoofVegPlant, dInt_dtRoofVegGround, IntRooftot, dInt_dtRooftot, dVRoofSoilVeg_dt, fRoofVeg, VRoofSoilVeg, OwRoofSoilVeg, OSwRoofSoilVeg, ExWaterRoofVeg_H, SoilPotWRoofVeg_H, SoilPotWRoofVeg_L, ExWaterRoofVeg_L, CiCO2LeafRoofVegSun, CiCO2LeafRoofVegShd, WBRoofVegInVeg, WBRoofVegInGround, WBRoofVegSoil, EBRoofImp, EBRoofVeg, Yroof, WBRoofImp, WBRoofVeg, WBRoofTot = eb_wb_roof(
                TR,
                TB,
                model.variables.temperature.tempvec,
                model.forcing.meteorological,
                model.variables.waterflux.Interception,
                model.variables.waterflux.ExWater,
                model.variables.waterflux.Vwater,
                model.variables.waterflux.Owater,
                model.variables.waterflux.SoilPotW,
                model.variables.waterflux.CiCO2Leaf,
                model.variables.waterflux.Runon,
                model.parameters.urbangeometry,
                model.parameters.surfacefractions.roof,
                model.parameters.soil.roof,
                model.parameters.optical.roof,
                model.parameters.thermal.roof,
                model.parameters.vegetation.roof,
                model.forcing.meteorological,
                model.forcing.anthropogenic,
                ParCalculation,
                BEM_on,
                RESPreCalc,
                rsRoofPreCalc,
            )

            SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t, albedo_canyon, LWRin_t, LWRout_t, LWRabs_t, LWREB_t, HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, HfluxGround, EfluxGroundImp, EfluxGroundBarePond, EfluxGroundBareSoil, EfluxGroundVegInt, EfluxGroundVegPond, EfluxGroundVegSoil, TEfluxGroundVeg, EfluxTreeInt, TEfluxTree, EfluxGroundBare, EfluxGroundVeg, EfluxGround, EfluxTree, LEfluxGroundImp, LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg, LEfluxTreeInt, LTEfluxTree, LEfluxGroundBare, LEfluxGroundVeg, LEfluxGround, LEfluxTree, CiCO2LeafTreeSun, CiCO2LeafTreeShd, CiCO2LeafGroundVegSun, CiCO2LeafGroundVegShd, raCanyontoAtm, raCanyontoAtmOrig, rap_can, rap_Htree_In, rb_HGround, rb_LGround, r_soilGroundbare, r_soilGroundveg, alp_soilGroundbare, alp_soilGroundveg, rs_sunGround, rs_shdGround, rs_sunTree, rs_shdTree, Fsun_L, Fshd_L, dw_L, RES_w1, RES_w2, rap_W1_In, rap_W2_In, rap_Zp1, HfluxWallSun, HfluxWallShade, EfluxWallSun, EfluxWallShade, LEfluxWallSun, LEfluxWallShade, HfluxCanyon, LEfluxCanyon, EfluxCanyon, G1WallSun, G2WallSun, dsWallSun, G1WallShade, G2WallShade, dsWallShade, G1GroundImp, TDampGroundImp, G1GroundBare, TDampGroundBare, G1GroundVeg, TDampGroundVeg, GTree, TDampTree, G1Ground, G1Canyon, G2Canyon, dsGroundImp, dsGroundBare, dsGroundVeg, dsTree, dsCanyonAir, Ycanyon, QTree, IntTree, dInt_dtTree, QGroundVegDrip, IntGroundVegPlant, dInt_dtGroundVegPlant, QGroundImp, IntGroundImp, dInt_dtGroundImp, fGroundImp, QGroundBarePond, IntGroundBare, dInt_dtGroundBare, fGroundBare, QGroundVegPond, IntGroundVegGround, dInt_dtGroundVegGround, fGroundVeg, VGroundSoilImp, OwGroundSoilImp, OSwGroundSoilImp, LkGroundImp, SoilPotWGroundImp_H, SoilPotWGroundImp_L, ExWaterGroundImp_H, ExWaterGroundImp_L, Rd_gimp, TEgveg_imp, TEtree_imp, Egimp_soil, dVGroundSoilImp_dt, Psi_Soil_gimp, Kf_gimp, VGroundSoilBare, OwGroundSoilBare, OSwGroundSoilBare, LkGroundBare, SoilPotWGroundBare_H, SoilPotWGroundBare_L, ExWaterGroundBare_H, ExWaterGroundBare_L, QGroundBareSoil, TEgveg_bare, TEtree_bare, Egbare_Soil, dVGroundSoilBare_dt, Psi_soil_gbare, Kf_gbare, VGroundSoilVeg, OwGroundSoilVeg, OSwGroundSoilVeg, LkGroundVeg, SoilPotWGroundVeg_H, SoilPotWGroundVeg_L, ExWaterGroundVeg_H, ExWaterGroundVeg_L, QGroundVegSoil, TEgveg_veg, TEtree_veg, Egveg_Soil, dVGroundSoilVeg_dt, Psi_soil_gveg, Kf_gveg, Qin_imp, Qin_bare, Qin_veg, Qin_bare2imp, Qin_bare2veg, Qin_imp2bare, Qin_imp2veg, Qin_veg2imp, Qin_veg2bare, VGroundSoilTot, OwGroundSoilTot, OSwGroundSoilTot, LkGround, Rd, dVGroundSoilTot_dt, SoilPotWGroundTot_L, ExWaterGroundTot_L, TEgveg_tot, SoilPotWGroundTot_H, ExWaterGroundTot_H, TEtree_tot, EB_TEtree, EB_TEgveg, WBIndv, WBTot, RunoffGroundTot, RunonGroundTot, Etot, DeepGLk, StorageTot, EBGroundImp, EBGroundBare, EBGroundVeg, EBTree, EBWallSun, EBWallShade, EBWallSunInt, EBWallShadeInt, EBCanyonT, EBCanyonQ, HumidityCan, HumidityAtm, u_Hcan, u_Zref_und, T2m, q2m, e_T2m, RH_T2m, qcan, e_Tcan, RH_Tcan, DHi, Himp_2m, Hbare_2m, Hveg_2m, Hwsun_2m, Hwshade_2m, Hcan_2m, DEi, Eimp_2m, Ebare_soil_2m, Eveg_int_2m, Eveg_soil_2m, TEveg_2m, Ecan_2m, dS_H_air, dS_LE_air = eb_wb_canyon(
                TC,
                TB,
                model.variables.temperature.tempvec,
                model.variables.humidity.Humidity,
                model.forcing.meteorological,
                model.variables.waterflux.Interception,
                model.variables.waterflux.ExWater,
                model.variables.waterflux.Vwater,
                model.variables.waterflux.Owater,
                model.variables.waterflux.SoilPotW,
                model.variables.waterflux.CiCO2Leaf,
                model.variables.temperature.tempdamp,
                model.variables.waterflux.Runon,
                model.variables.waterflux.Qinlat,
                ViewFactor,
                model.parameters.urbangeometry,
                model.parameters.surfacefractions.ground,
                WallLayers,
                model.parameters.soil.ground,
                ParInterceptionTree,
                model.parameters.optical.ground,
                model.parameters.optical.wall,
                model.parameters.optical.tree,
                model.parameters.thermal.ground,
                model.parameters.thermal.wall,
                model.parameters.vegetation.ground,
                model.parameters.vegetation.tree,
                model.forcing.sunposition,
                model.forcing.meteorological,
                model.forcing.anthropogenic,
                ParCalculation,
                model.variables.buildingenergymodel.TempVecB,
                G2Roof,
                model.parameters.building_energy.indoor_optical,
                ParHVAC,
                model.parameters.building_energy.thermal,
                model.parameters.building_energy.windows,
                BEM_on,
                RESPreCalc,
                fconvPreCalc,
                fconv,
                rsGroundPreCalc,
                rsTreePreCalc,
                model.forcing.hvacschedule,
            )

            SWRinWsun = SWRabs_t.WallSun
            SWRinWshd = SWRabs_t.WallShade

            HbuildIntc, LEbuildIntc, GbuildIntc, SWRabsBc, LWRabsBc, TDampGroundBuild, WasteHeat, EnergyUse, HumidBuilding, ParACHeat_t, YBuildInt = BuildingEnergyModel.eb_solver_building_output(
                TC,
                TB,
                model.variables.buildingenergymodel.TempVecB,
                model.variables.temperature.tempvec,
                model.variables.humidity.Humidity,
                model.forcing.meteorological,
                SWRinWsun,
                SWRinWshd,
                G2Roof,
                G2WallSun,
                G2WallShade,
                model.variables.temperature.tempdamp,
                SWRabs_t,
                model.parameters.urbangeometry,
                model.parameters.building_energy.indoor_optical,
                ParHVAC,
                ParCalculation,
                model.parameters.building_energy.thermal,
                model.parameters.building_energy.windows,
                BEM_on,
                model.forcing.hvacschedule,
            )

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
            @infiltrate
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
            # Need to update all fields
            ExWater_t.VRoofSoilVeg[:] = VRoofSoilVeg
            ExWater_t.VGroundSoilImp[:] = VGroundSoilImp
            ExWater_t.VGroundSoilBare[:] = VGroundSoilBare
            ExWater_t.VGroundSoilVeg[:] = VGroundSoilVeg
            ExWater_t.VGroundSoilTot[:] = VGroundSoilTot

            # Owater
            # Owater_t was never created?
            Owater_t.OwRoofSoilVeg[:] = OwRoofSoilVeg
            Owater_t.OwGroundSoilImp[:] = OwGroundSoilImp
            Owater_t.OwGroundSoilBare[:] = OwGroundSoilBare
            Owater_t.OwGroundSoilVeg[:] = OwGroundSoilVeg
            OWater_t.OwGroundSoilTot[:] = OwGroundSoilTot

            # SoilPotW
            SoilPotW_t.SoilPotWRoofVeg_H[:] = SoilPotWRoofVeg_H[] # technically a matrix
            SoilPotW_t.SoilPotWRoofVeg_L[:] = SoilPotWRoofVeg_L[] # vector
            SoilPotW_t.SoilPotWGroundImp_H[:] = SoilPotWGroundImp_H
            SoilPotW_t.SoilPotWGroundImp_L[:] = SoilPotWGroundImp_L
            SoilPotW_t.SoilPotWGroundBare_H[:] = SoilPotWGroundBare_H
            SoilPotW_t.SoilPotWGroundBare_L[:] = SoilPotWGroundBare_L
            SoilPotW_t.SoilPotWGroundVeg_H[:] = SoilPotWGroundVeg_H
            SoilPotW_t.SoilPotWGroundVeg_L[:] = SoilPotWGroundVeg_L
            SoilPotW_t.SoilPotWGroundTot_L[:] = SoilPotWGroundTot_L
            SoilPotW_t.SoilPotWGroundTot_H[:] = SoilPotWGroundTot_H

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
            Runon_t.RunonUrbanTot = urban_average(
                RunonRoofTot, RunonGroundTot, model.parameters.urbangeometry
            )
            Runon_t.RunoffUrbanTot = urban_average(
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
        end

        Tmrt, BoleanInSun, SWRdir_Person, SWRdir_in_top, SWRdir_in_bottom, SWRdir_in_east, SWRdir_in_south, SWRdir_in_west, SWRdir_in_north, SWRdiff_Person, LWR_Person = MeanRadiantTemperature.mean_radiant_temperature(
            SWRout_t,
            LWRout_t,
            model.forcing.meteorological,
            ViewFactorPoint,
            model.parameters.vegetation.tree,
            model.parameters.urbangeometry,
            model.forcing.sunposition,
            model.parameters.person,
            FT(Dates.hour(model.forcing.datetime[])),
        )
        # TODO: check whether we should be using the hour as a float (e.g. 10.5 for 10:30) or
        # the hour as an integer (10 for 10:30).

        u_ZPerson = Resistance.wind_profile_point_output(
            model.parameters.person.HeightWind,
            model.parameters.urbangeometry,
            model.parameters.vegetation.tree,
            model.forcing.meteorological,
            model.parameters.surfacefractions.ground,
            model.parameters.vegetation.ground,
        )

        UTCI_approx = OutdoorThermalComfort.utci_approx(
            T2m - FT(273.15), RH_T2m * 100, Tmrt, u_ZPerson
        )

        # Assign outputs
        # Urban average store in *TotalUrban field, which is not part of the RadiationFluxes
        # composite type

        # tempvec - already done higher

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
        update!(model.variables.waterflux.Vwater, ExWater_t)
        update!(
            model.variables.waterflux.Owater,
            Owater_t,
            model.parameters.soil.roof,
            model.parameters.soil.ground,
        )
        update!(model.variables.waterflux.SoilPotW, SoilPotW_t)
        update!(model.variables.waterflux.CiCO2Leaf, CiCO2Leaf_t)
        update!(model.variables.waterflux.Runon, Runon_t)
        update!(model.variables.waterflux.Qinlat, Qinlat_t)
    end
    return Ttot, fval, exitflag
end

function update!(
    x::ModelComponents.ModelVariables.TempVec{FT}, Ttot::Vector{FT}
) where {FT<:AbstractFloat}
    x.TRoofImp = Ttot[1]
    x.TRoofVeg = Ttot[2]
    x.TRoofIntImp = Ttot[3]
    x.TRoofIntVeg = Ttot[4]
    x.TGroundImp = Ttot[5]
    x.TGroundBare = Ttot[6]
    x.TGroundVeg = Ttot[7]
    x.TWallSun = Ttot[8]
    x.TWallShade = Ttot[9]
    x.TTree = Ttot[10]
    x.TWallIntSun = Ttot[11]
    x.TWallIntShade = Ttot[12]
    x.TCanyon = Ttot[13]

    return nothing
end

function update!(
    x::ModelComponents.ModelVariables.Humidity{FT}, Humiditytot::Vector{FT}
) where {FT<:AbstractFloat}
    x.CanyonSpecific = Humiditytot[14]
    return nothing
end

function update!(
    x::ModelComponents.ModelVariables.TempVecB{FT}, Ttot::Vector{FT}
) where {FT<:AbstractFloat}
    x.Tceiling = Ttot[15]
    x.Tinwallsun = Ttot[16]
    x.Tinwallshd = Ttot[17]
    x.Twindows = Ttot[18]
    x.Tinground = Ttot[19]
    x.Tintmass = Ttot[20]
    x.Tbin = Ttot[21]
    x.qbin = Ttot[22]
    return nothing
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
