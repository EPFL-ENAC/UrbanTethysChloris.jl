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
        dts=3600.0,
        row=1000.0,
        cp_atm=model.forcing.meteorological.cp_atm,
        rho_atm=model.forcing.meteorological.rho_atm,
    )

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

        for HVACittm in 1:2
            if BEM_on && HVACittm == 2
                if !ParHVACorig.Acon && !ParHVACorig.Heatingon
                    continue
                end

                ## Add complete logic with EnergyUse (named tuple created at end of first iteration)
            end
            @info RESPreCalc fconvPreCalc BEM_on fconv

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
        end
    end
end
