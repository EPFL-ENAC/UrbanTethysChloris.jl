using UrbanTethysChloris
using YAML
using NCDatasets

FT = Float32

yaml_path = joinpath(@__DIR__, "data", "parameters.yaml")
ncdf_path = joinpath(@__DIR__, "data", "input_data.nc")

model, forcing = create_model(FT, ncdf_path, yaml_path);

initialize!(model, forcing)
