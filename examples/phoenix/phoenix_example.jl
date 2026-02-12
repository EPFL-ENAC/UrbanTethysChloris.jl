using UrbanTethysChloris
using YAML
using NCDatasets

FT = Float64

yaml_path = joinpath(@__DIR__, "data", "phoenix_parameters.yaml")
ncdf_path = joinpath(@__DIR__, "data", "phoenix_data.nc")

model, forcing = create_model(FT, ncdf_path, yaml_path);

initialize!(model, forcing)

O33 = (
    roof=model.variables.waterflux.Owater.OwRoofSoilVeg[1],
    ground=model.variables.waterflux.Owater.OwGroundSoilVeg[1],
)

NN = 100
results, view_factor_out, view_factor_point_out = run_simulation(
    model,
    forcing;
    NN=NN,
    O33=O33,
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
