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

model, forcing = create_model(FT, ncdf_path, yaml_path);
O33 = (
    roof=model.variables.waterflux.Owater.OwRoofSoilVeg[1],
    ground=model.variables.waterflux.Owater.OwGroundSoilVeg[1],
)

initialize!(model, forcing)

NN = 200
results, view_factor_out, view_factor_point_out = run_simulation(
    model,
    forcing;
    NN=NN,
    O33=O33,
    ViewFactors=(view_factor, view_factor_point),
    fconvPreCalc=false,
    output_level=extended_outputs,
    WallLayers=(dz1_wall=0.11, dz2_wall=0.11),
)

x, x_day, x_month, fig1, fig2 = urban_climate_variables(results, model, forcing, NN)

ef_urban, ef_canyon, ef_roof, fig3, fig4, fig5 = plan_area_energy_balance_calculation(
    results, model, forcing, view_factor_out, NN
)

wf_urban, wf_canyon, wf_roof, wf_building, fig6, fig7 = water_balance_components(
    results, model, forcing, NN
)
