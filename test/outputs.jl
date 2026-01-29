using Test
using UrbanTethysChloris.ModelComponents: WaterFluxVariables, ExtendedOutputs, accessors
using UrbanTethysChloris.ModelComponents.Parameters: initialize_soil_parameters
using ..TestUtils: load_test_parameters

FT = Float64

input_data = load_test_parameters()

soil = initialize_soil_parameters(Float64, input_data["soil"])

MR = soil.roof.ms
MG = soil.ground.ms
water_flux_vars = WaterFluxVariables(FT, soil)

x = accessors(WaterFluxVariables, ExtendedOutputs)

# How the accessor works in practice
water_flux_vars.Runoff.QTree = 17
@test x[:Runoff][:QTree](water_flux_vars) == 17
