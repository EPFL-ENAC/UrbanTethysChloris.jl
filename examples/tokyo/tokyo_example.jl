using UrbanTethysChloris
using YAML
using NCDatasets

FT = Float64

yaml_path = joinpath(@__DIR__, "data", "tokyo_parameters.yaml")
ncdf_path = joinpath(@__DIR__, "data", "tokyo_data.nc")

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

using MAT
using Plots

matdata = matread("testrun_tokyo_1000steps_20260203_1217.mat")

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

plot_differences(
    matdata["Humidity"]["CanyonSpecific"], results[:Humidity][:CanyonSpecific], 1:Nx
)

plot_differences(matdata["TempVecB"]["Tceiling"], results[:TempVecB][:Tceiling], 1:NN)

plot_differences(
    matdata["dInt_dt"]["dInt_dtGroundVegPlant"],
    results[:dInt_dt][:dInt_dtGroundVegPlant],
    1:NN,
)

plot_differences(
    matdata["MeanRadiantTemperature"]["Tmrt"], results[:mrt][:Tmrt], 1:NN; relative=false
)
