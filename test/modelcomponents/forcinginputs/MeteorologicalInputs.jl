using Test
using UrbanTethysChloris.ModelComponents.ForcingInputs:
    MeteorologicalInputs, initialize_meteorological_inputs
using ....TestUtils: load_test_netcdf

FT = Float64

ds = load_test_netcdf();
theta_Z = zeros(FT, ds.dim["hours"])

@test_nowarn initialize_meteorological_inputs(FT, ds, theta_Z);
