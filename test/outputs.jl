using Test
using UrbanTethysChloris.ModelComponents: WaterFluxVariables, ExtendedOutputs, accessors
using UrbanTethysChloris.ModelComponents.Parameters: initialize_soil_parameters
using UrbanTethysChloris: _create_model
using ..TestUtils: load_test_parameters, load_test_netcdf

FT = Float64

yaml_data = load_test_parameters();
netcdf_data = load_test_netcdf();

x = accessors(WaterFluxVariables, ExtendedOutputs)

model, forcing = _create_model(FT, netcdf_data, yaml_data);

# # How the accessor works in practice
model.variables.waterflux.Runoff.QTree = 17
@test x[:Runoff][:QTree](model) == 17
