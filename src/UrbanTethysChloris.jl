module UrbanTethysChloris

using TethysChlorisCore
using Dates
using Statistics
using ConstructionBase
using LeastSquaresOptim

include(joinpath("modelcomponents", "ModelComponents.jl"))
using .ModelComponents

include(joinpath("ray_tracing", "RayTracing.jl"))
using .RayTracing

include(joinpath("mean_radiant_temperature", "MeanRadiantTemperature.jl"))
using .MeanRadiantTemperature

include(joinpath("outdoor_thermal_comfort", "OutdoorThermalComfort.jl"))
using .OutdoorThermalComfort

include(joinpath("radiation", "Radiation.jl"))
using .Radiation

include(joinpath("soil", "Soil.jl"))
using .Soil

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

include("eb_solver_canyon.jl")
include("eb_solver_roof.jl")
include("eb_wb_roof.jl")
include("eb_solver_urban_climate_building_energy_model.jl")
include("f_solver_tot.jl")

end
