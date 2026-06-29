using UrbanTethysChloris
using YAML
using NCDatasets

FT = Float64

using UrbanTethysChloris.RayTracing: ViewFactor, ViewFactorPoint

view_factor = ViewFactor{FT}(;
    F_gs_nT=0.381290805156702,
    F_gw_nT=0.309354597421649,
    F_ww_nT=0.447921026139827,
    F_wg_nT=0.276039486930087,
    F_ws_nT=0.276039486930087,
    F_sg_nT=0.381290805156702,
    F_sw_nT=0.309354597421649,
    F_gs_T=0.290910447761194,
    F_gt_T=0.102572139303483,
    F_gw_T=0.303258706467662,
    F_ww_T=0.346598392652124,
    F_wt_T=0.147054726368159,
    F_wg_T=0.270600076540375,
    F_ws_T=0.235746804439342,
    F_sg_T=0.290910447761194,
    F_sw_T=0.264199004975124,
    F_st_T=0.180691542288557,
    F_tg_T=0.163248629936597,
    F_tw_T=0.262291660245468,
    F_ts_T=0.287579521301221,
    F_tt_T=0.024588528271246,
);

view_factor_point = ViewFactorPoint{Float64}(;
    F_pg=0.392328358208955,
    F_pwLeft=0.229149253731343,
    F_pwRight=0.221840796019900,
    F_ps=0.110054726368159,
    F_pt=0.046626865671642,
);

yaml_path = joinpath(@__DIR__, "data", "newdelhi_parameters.yaml")
ncdf_path = joinpath(@__DIR__, "data", "newdelhi_data.nc")

options = ModelOptions(; fconvPreCalc=false, output_level=extended_outputs)

model, forcing = create_model(FT, ncdf_path, yaml_path);

initialize!(model, forcing)

O33 = (
    roof=model.variables.waterflux.Owater.OwRoofSoilVeg[1],
    ground=model.variables.waterflux.Owater.OwGroundSoilVeg[1],
)

NN = 100
results, view_factor_out, view_factor_point_out = run_simulation(
    model, forcing; options=options, NN=NN, O33=O33
)

x, x_day, x_month, fig1, fig2 = urban_climate_variables(results, model, forcing, NN)

ef_urban, ef_canyon, ef_roof, fig3, fig4, fig5 = plan_area_energy_balance_calculation(
    results, model, forcing, view_factor_out, NN
)

wf_urban, wf_canyon, wf_roof, wf_building, fig6, fig7 = water_balance_components(
    results, model, forcing, NN
)

# Simulation starts a bit different from MATLAB from the first iteration
using MAT
using Plots

matdata = matread("testrun_newdelhi_1000steps_20260203_1046.mat")

function plot_differences(matlab_variable, julia_variable, NNv; relative=false)
    p1 = plot(
        NNv,
        [matlab_variable[NNv], julia_variable[NNv]];
        label=["MATLAB" "Julia"],
        color=[:blue :red],
        grid=true,
    )

    err = julia_variable[NNv] - matlab_variable[NNv]
    err_label = "Difference"
    if relative
        err ./= matlab_variable[NNv]
        err_label = "Relative Difference"
    end

    p2 = plot(NNv, err; label=err_label, color=:green, grid=true)

    return plot(p1, p2; layout=(2, 1))
end

Nx = NN
# Somewhat big differences in temperature results
plot_differences(matdata["TempVec"]["TCanyon"], results[:tempvec][:TCanyon], 1:Nx)

plot_differences(matdata["TempVec"]["TRoofImp"], results[:tempvec][:TRoofImp], 1:Nx)

plot_differences(matdata["TempVecB"]["Tceiling"], results[:TempVecB][:Tceiling], 1:NN)

plot_differences(
    matdata["dInt_dt"]["dInt_dtGroundVegPlant"],
    results[:dInt_dt][:dInt_dtGroundVegPlant],
    1:Nx,
)

plot_differences(
    matdata["MeanRadiantTemperature"]["Tmrt"], results[:mrt][:Tmrt], 1:Nx; relative=false
)

plot_differences(matdata["EB"]["EBCanyonQ"], results[:EB][:EBCanyonQ], 1:Nx)
# Instability comes from mostly from WallSun and WallShade

# SWRabs is fine
plot_differences(matdata["SWRabs"]["SWRabsWallSunExt"], results[:SWRabs][:WallSunExt], 1:47)

# The difference appears to come from these three variables
plot_differences(matdata["LWRabs"]["LWRabsWallSun"], results[:LWRabs][:WallSun], 1:47)
# Difference reaches 0.08 at iteration 46

plot_differences(matdata["Gflux"]["G1WallSun"], results[:Gflux][:G1WallSun], 1:47)
# Difference reaches 0.06 at iteration 46

plot_differences(matdata["Hflux"]["HfluxWallSun"], results[:Hflux][:HfluxWallSun], 1:47)
# Difference reaches 0.20 at iteration 46
# The difference most likely comes from the temperature itself, with differences reaching
# up to 0.05
# ExWater variables are quite (very!) different, despite their amplitude
# Vwater, Owater, SoilPotW are different by 0.03 for the ground, not for the roof
# These values are normal given the low convergence threshold
plot_differences(
    matdata["ExWater"]["ExWaterGroundTot_L"], results[:ExWater][:ExWaterGroundTot_L], 1:Nx
)

plot_differences(matdata["EB"]["EBRoofVeg"], results[:EB][:EBRoofVeg], 1:Nx; relative=true)

# There's a weird pattern with the initial value of EBRoofVeg being extremely high at the first timestep
