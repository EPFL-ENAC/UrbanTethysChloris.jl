using UrbanTethysChloris
using YAML
using NCDatasets

yaml_path = joinpath(@__DIR__, "data", "parameters.yaml")
ncdf_path = joinpath(@__DIR__, "data", "input_data.nc")

model = initialize_model(Float32, ncdf_path, yaml_path);
