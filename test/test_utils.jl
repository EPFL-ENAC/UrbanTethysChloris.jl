module TestUtils

using NCDatasets
using YAML

"""
    load_test_netcdf(path::String = joinpath(@__DIR__, "data", "input_data.nc"))

Load NetCDF input data for testing.

# Arguments
- `path`: Path to NetCDF file containing test data, defaults to data/input_data.nc

# Returns
- `NCDataset`: Dataset containing test input data
"""
function load_test_netcdf(path::String=joinpath(@__DIR__, "data", "input_data.nc"))
    return NCDataset(path)
end

"""
    load_test_parameters(path::String = joinpath(@__DIR__, "data", "parameters.yaml"))

Load and preprocess model parameters from YAML file for testing.

# Arguments
- `path`: Path to YAML parameter file, defaults to data/parameters.yaml

# Returns
- `Dict`: Preprocessed parameter dictionary with:
"""
function load_test_parameters(path::String=joinpath(@__DIR__, "data", "parameters.yaml"))
    return YAML.load_file(path; dicttype=Dict{String,Any})
end

end
