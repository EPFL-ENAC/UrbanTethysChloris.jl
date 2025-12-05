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
    ViewFactor::RayTracing.ViewFactor{FT}=nothing,
) where {FT<:AbstractFloat}

    # view factor
    if isnothing(ViewFactor)
        ViewFactor, ViewFactorPoint = RayTracing.view_factors_canyon(
            model.parameters.urbangeometry, model.parameters.person, mc_sample_size, n_rays
        )
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

    Ttot, fval, exitflag = nothing, nothing, nothing

    for i in 1:NN
        @info "Starting iteration $i / $NN"

        # Hard-code ittm2Ext with values from first iteration, keep as named tuple for dev
        TempVec_ittm2Ext = (;
            TRoofImp=fill(model.variables.temperature.tempvec.TRoofImp, 4),
            TRoofVeg=fill(model.variables.temperature.tempvec.TRoofVeg, 4),
            TRoofIntImp=fill(model.variables.temperature.tempvec.TRoofIntImp, 4),
            TRoofIntVeg=fill(model.variables.temperature.tempvec.TRoofIntVeg, 4),
            TGroundImp=fill(model.variables.temperature.tempvec.TGroundImp, 4),
            TGroundBare=fill(model.variables.temperature.tempvec.TGroundBare, 4),
            TGroundVeg=fill(model.variables.temperature.tempvec.TGroundVeg, 4),
            TTree=fill(model.variables.temperature.tempvec.TTree, 4),
            TWallSun=fill(model.variables.temperature.tempvec.TWallSun, 4),
            TWallShade=fill(model.variables.temperature.tempvec.TWallShade, 4),
            TWallIntSun=fill(model.variables.temperature.tempvec.TWallIntSun, 4),
            TWallIntShade=fill(model.variables.temperature.tempvec.TWallIntShade, 4),
            TCanyon=fill(model.variables.temperature.tempvec.TCanyon, 4),
            Tatm=fill(model.variables.temperature.tempvec.Tatm, 4),
        )

        Humidity_ittm2Ext = (;
            CanyonRelative=fill(model.variables.humidity.Humidity.CanyonRelative, 4),
            CanyonSpecific=fill(model.variables.humidity.Humidity.CanyonSpecific, 4),
            CanyonVapourPre=fill(model.variables.humidity.Humidity.CanyonVapourPre, 4),
            CanyonRelativeSat=fill(model.variables.humidity.Humidity.CanyonRelativeSat, 4),
            CanyonSpecificSat=fill(model.variables.humidity.Humidity.CanyonSpecificSat, 4),
            CanyonVapourPreSat=fill(
                model.variables.humidity.Humidity.CanyonVapourPreSat, 4
            ),
            AtmRelative=fill(model.variables.humidity.Humidity.AtmRelative, 4),
            AtmSpecific=fill(model.variables.humidity.Humidity.AtmSpecific, 4),
            AtmVapourPre=fill(model.variables.humidity.Humidity.AtmVapourPre, 4),
            AtmRelativeSat=fill(model.variables.humidity.Humidity.AtmRelativeSat, 4),
            AtmSpecificSat=fill(model.variables.humidity.Humidity.AtmSpecificSat, 4),
            AtmVapourPreSat=fill(model.variables.humidity.Humidity.AtmVapourPreSat, 4),
        )

        TempVecB_ittm2Ext = (;
            Tceiling=fill(model.variables.buildingenergymodel.TempVecB.Tceiling, 4),
            Tinwallsun=fill(model.variables.buildingenergymodel.TempVecB.Tinwallsun, 4),
            Tinwallshd=fill(model.variables.buildingenergymodel.TempVecB.Tinwallshd, 4),
            Twindows=fill(model.variables.buildingenergymodel.TempVecB.Twindows, 4),
            Tinground=fill(model.variables.buildingenergymodel.TempVecB.Tinground, 4),
            Tintmass=fill(model.variables.buildingenergymodel.TempVecB.Tintmass, 4),
            Tbin=fill(model.variables.buildingenergymodel.TempVecB.Tbin, 4),
            qbin=fill(model.variables.buildingenergymodel.TempVecB.qbin, 4),
        )

        Meteo_ittm = (; SWRin=zeros(FT, 2), Rain=zeros(FT, 2))

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
                @infiltrate
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

            SWRabsRoofImp, SWRabsRoofVeg, SWRabsTotalRoof, SWRoutRoofImp, SWRoutRoofVeg, SWRoutTotalRoof, SWRinRoofImp, SWRinRoofVeg, SWRinTotalRoof, SWREBRoofImp, SWREBRoofVeg, SWREBTotalRoof, LWRabsRoofVeg, LWRabsRoofImp, LWRabsTotalRoof, LWRoutRoofVeg, LWRoutRoofImp, LWRoutTotalRoof, LWRinRoofImp, LWRinRoofVeg, LWRinTotalRoof, LWREBRoofImp, LWREBRoofVeg, LWREBTotalRoof, HfluxRoofImp, HfluxRoofVeg, HfluxRoof, LEfluxRoofImp, LEfluxRoofVegInt, LEfluxRoofVegPond, LEfluxRoofVegSoil, LTEfluxRoofVeg, LEfluxRoofVeg, LEfluxRoof, G1RoofImp, G2RoofImp, dsRoofImp, G1RoofVeg, G2RoofVeg, dsRoofVeg, G1Roof, G2Roof, dsRoof, raRooftoAtm, rb_LRoof, rap_LRoof, r_soilRoof, rs_sunRoof, rs_shdRoof, EfluxRoofImp, EfluxRoofVegInt, EfluxRoofVegPond, EfluxRoofVegSoil, TEfluxRoofVeg, EfluxRoofVeg, EfluxRoof, QRoofImp, QRoofVegDrip, QRoofVegPond, LkRoofImp, LkRoofVeg, LkRoof, QRoofVegSoil, RunoffRoofTot, RunonRoofTot, IntRoofImp, IntRoofVegPlant, IntRoofVegGround, dInt_dtRoofImp, dInt_dtRoofVegPlant, dInt_dtRoofVegGround, IntRooftot, dInt_dtRooftot, dVRoofSoil_dt, fRoofVeg, VRoofSoil, OwRoofSoil, OSwRoofSoil, ExWaterRoof_H, SoilPotWRoof_H, ExWaterRoof_L, SoilPotWRoof_L, CiCO2LeafRoofVegSun, CiCO2LeafRoofVegShd, WBRoofVegInVeg, WBRoofVegInGround, WBRoofVegSoil, EBRoofImp, EBRoofVeg, Yroof, WBRoofImp, WBRoofVeg, WBRoofTot = eb_wb_roof(
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

            SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t, albedo_canyon, LWRin_t, LWRout_t, LWRabs_t, LWREB_t, HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, HfluxGround, EfluxGroundImp, EfluxGroundBarePond, EfluxGroundBareSoil, EfluxGroundVegInt, EfluxGroundVegPond, EfluxGroundVegSoil, TEfluxGroundVeg, EfluxTreeInt, TEfluxTree, EfluxGroundBare, EfluxGroundVeg, EfluxGround, EfluxTree, LEfluxGroundImp, LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg, LEfluxTreeInt, LTEfluxTree, LEfluxGroundBare, LEfluxGroundVeg, LEfluxGround, LEfluxTree, CiCO2LeafTreeSun, CiCO2LeafTreeShd, CiCO2LeafGroundVegSun, CiCO2LeafGroundVegShd, raCanyontoAtm, raCanyontoAtmOrig, rap_can, rap_Htree_In, rb_HGround, rb_LGround, r_soilGroundbare, r_soilGroundveg, alp_soil_bare, alp_soil_veg, rs_sunGround, rs_shdGround, rs_sunTree, rs_shdTree, Fsun_L, Fshd_L, dw_L, RES_w1, RES_w2, rap_W1_In, rap_W2_In, rap_Zp1, HfluxWallSun, HfluxWallShade, EfluxWallSun, EfluxWallShade, LEfluxWallSun, LEfluxWallShade, HfluxCanyon, LEfluxCanyon, EfluxCanyon, G1WallSun, G2WallSun, dsWallSun, G1WallShade, G2WallShade, dsWallShade, G1GroundImp, TDampGroundImp, G1GroundBare, TDampGroundBare, G1GroundVeg, TDampGroundVeg, GTree, TDampTree, G1Ground, G1Canyon, G2Canyon, dsGroundImp, dsGroundBare, dsGroundVeg, dsTree, dsCanyonAir, Ycanyon, q_tree_dwn, In_tree, dIn_tree_dt, q_gveg_dwn, In_gveg, dIn_gveg_dt, q_gimp_runoff, In_gimp, dIn_gimp_dt, f_inf_gimp, q_gbare_runoff, In_gbare, dIn_gbare_dt, f_inf_gbare, q_gveg_runoff, In_gveg_pond, dIn_gveg_pond_dt, f_inf_gveg, V_gimp, O_gimp, OS_gimp, Lk_gimp, Psi_s_H_gimp, Psi_s_L_gimp, Exwat_H_gimp, Exwat_L_gimp, Rd_gimp, TEgveg_imp, TEtree_imp, Egimp_soil, dV_dt_gimp, Psi_Soil_gimp, Kf_gimp, V_gbare, O_gbare, OS_gbare, Lk_gbare, Psi_s_H_gbare, Psi_s_L_gbare, Exwat_H_gbare, Exwat_L_gbare, Rd_gbare, TEgveg_bare, TEtree_bare, Egbare_Soil, dV_dt_gbare, Psi_soil_gbare, Kf_gbare, V_gveg, O_gveg, OS_gveg, Lk_gveg, Psi_s_H_gveg, Psi_s_L_gveg, Exwat_H_gveg, Exwat_L_gveg, Rd_gveg, TEgveg_veg, TEtree_veg, Egveg_Soil, dV_dt_gveg, Psi_soil_gveg, Kf_gveg, Qin_imp, Qin_bare, Qin_veg, Qin_bare2imp, Qin_bare2veg, Qin_imp2bare, Qin_imp2veg, Qin_veg2imp, Qin_veg2bare, V, O, OS, Lk, Rd, dV_dt, Psi_s_L, Exwat_L, TEgveg_tot, Psi_s_H_tot, Exwat_H, TEtree_tot, EB_TEtree, EB_TEgveg, WBIndv, WBTot, Runoff, Runon_ittm, Etot, DeepGLk, StorageTot, EBGroundImp, EBGroundBare, EBGroundVeg, EBTree, EBWallSun, EBWallShade, EBWallSunInt, EBWallShadeInt, EBCanyonT, EBCanyonQ, HumidityCan, HumidityAtm, u_Hcan, u_Zref_und, T2m, q2m, e_T2m, RH_T2m, qcan, e_Tcan, RH_Tcan, DHi, Himp_2m, Hbare_2m, Hveg_2m, Hwsun_2m, Hwshade_2m, Hcan_2m, DEi, Eimp_2m, Ebare_soil_2m, Eveg_int_2m, Eveg_soil_2m, TEveg_2m, Ecan_2m, dS_H_air, dS_LE_air = eb_wb_canyon(
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

            HbuildInt, LEbuildInt, GbuildInt, SWRabsB, LWRabsB, Tdpfloor, WasteHeat, EnergyUse, HumidityBuilding, ParACHeat, YBuildInt = BuildingEnergyModel.eb_solver_building_output(
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
        end

        # MRT

        # Wind profile output

        # UTCI

        # Assign outputs
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
