using UrbanTethysChloris
using UrbanTethysChloris.RayTracing: ViewFactor, ViewFactorPoint
using YAML
using NCDatasets

FT = Float64

view_factor = ViewFactor{FT}(;
    F_gs_nT=0.38030860180886683,
    F_gw_nT=0.3098456990955666,
    F_ww_nT=0.4489513413008077,
    F_wg_nT=0.27552432934959614,
    F_ws_nT=0.27552432934959614,
    F_sg_nT=0.38030860180886683,
    F_sw_nT=0.3098456990955666,
    F_gs_T=0.3203681592039802,
    F_gt_T=0.09739303482587065,
    F_gw_T=0.29111940298507455,
    F_ww_T=0.3613650593187906,
    F_wt_T=0.1503432835820897,
    F_wg_T=0.25887233065442017,
    F_ws_T=0.2294193264446995,
    F_sg_T=0.3203681592039802,
    F_sw_T=0.2579975124378108,
    F_st_T=0.16363681592039817,
    F_tg_T=0.15500582915258426,
    F_tw_T=0.2690851190794379,
    F_ts_T=0.2604360812554992,
    F_tt_T=0.04638785143304075,
);

view_factor_point = ViewFactorPoint{Float64}(;
    F_pg=0.38348258706467636,
    F_pwLeft=0.21121890547263672,
    F_pwRight=0.2085671641791044,
    F_ps=0.15139303482587063,
    F_pt=0.045338308457711474,
);

yaml_path = joinpath(@__DIR__, "data", "parameters.yaml")
ncdf_path = joinpath(@__DIR__, "data", "input_data.nc")

model, forcing = create_model(FT, ncdf_path, yaml_path);
O33 = (
    roof=model.variables.waterflux.Owater.OwRoofSoilVeg[1],
    ground=model.variables.waterflux.Owater.OwGroundSoilVeg[1],
)

initialize!(model, forcing)

NN = 100
results = run_simulation(
    model,
    forcing;
    NN=NN,
    ViewFactors=(view_factor, view_factor_point),
    O33=O33,
    fconvPreCalc=false,
)

x, x_day, x_month, fig1, fig2 = urban_climate_variables(results, model, forcing, NN)
ef_urban, ef_canyon, ef_roof, fig3, fig4, fig5 = plan_area_energy_balance_calculation(
    results, model, forcing, NN
)

wf_urban, wf_canyon, wf_roof, wf_building, fig6, fig7 = water_balance_components(
    results, model, forcing, NN
)
