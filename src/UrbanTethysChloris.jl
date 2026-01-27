module UrbanTethysChloris

using TethysChlorisCore
using Dates
using Statistics
using ConstructionBase
using LeastSquaresOptim
using Roots
using NCDatasets
using YAML
using NaNMath: NaNMath
using StaticArraysCore
using DataFrames: DataFrame, select, combine, groupby, Not
using DataFramesMeta: @chain
using Plots: plot, plot!
using LaTeXStrings: @L_str

include(joinpath("soil", "Soil.jl"))
using .Soil

include(joinpath("modelcomponents", "ModelComponents.jl"))
using .ModelComponents
import .ModelComponents: update!

include("Model.jl")
export create_model, initialize!

include(joinpath("ray_tracing", "RayTracing.jl"))
using .RayTracing

include(joinpath("outdoor_thermal_comfort", "OutdoorThermalComfort.jl"))
using .OutdoorThermalComfort

include(joinpath("radiation", "Radiation.jl"))
using .Radiation

# Depends on Radiation, RayTracing
include(joinpath("mean_radiant_temperature", "MeanRadiantTemperature.jl"))
using .MeanRadiantTemperature

include(joinpath("conductive_heat", "ConductiveHeat.jl"))
using .ConductiveHeat

include(joinpath("building_energy_model", "BuildingEnergyModel.jl"))
using .BuildingEnergyModel

include(joinpath("resistance", "Resistance.jl"))
using .Resistance

include(joinpath("turbulent_heat", "TurbulentHeat.jl"))
using .TurbulentHeat

include(joinpath("water", "Water.jl"))
using .Water

include("incoming_longwave.jl")
include("set_sun_variables.jl")

include("extrapolate.jl")
export update!

include("eb_solver_canyon.jl")
include("eb_solver_roof.jl")
include("eb_wb_roof.jl")
include("eb_solver_urban_climate_building_energy_model.jl")
include("f_solver_tot.jl")
include("eb_wb_canyon.jl")
include("update_hvac_parameters.jl")
include("run_simulation.jl")
export run_simulation

include("urban_climate_variables.jl")
export urban_climate_variables

include("plan_area_energy_balance_calculation.jl")
export plan_area_energy_balance_calculation

include("post_calculate_soil_moisture_change.jl")
include("water_balance_components.jl")
export water_balance_components

end
