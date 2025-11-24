using UrbanTethysChloris
using YAML
using NCDatasets

yaml_path = joinpath(@__DIR__, "data", "parameters.yaml")
ncdf_path = joinpath(@__DIR__, "data", "input_data.nc")

model, forcing = create_model(Float32, ncdf_path, yaml_path);
