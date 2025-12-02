
function run_simulation(
    model::Model;
    NN::Signed=nothing,
    mc_sample_size::Int=1000,
    n_rays::Int=200,
    RESPreCalc::Bool=true,
    fconvPreCalc::Bool=true,
    BEM_on::Bool=true,
)

    # view factor
    ViewFactor, ViewFactorPoint = RayTracing.view_factors_canyon(
        model.parameters.urbangeometry, model.parameters.person, mc_sample_size, n_rays
    )

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

    for i in 1:NN
        @info "Starting iteration $i / $NN"
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
    end
end
