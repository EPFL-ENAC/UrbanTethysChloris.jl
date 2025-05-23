module Resistance

using Roots
using ..Soil

include("aerodynamic_resistance.jl")
include("backcalculate_obhukov_length.jl")
include("canopy_resistance_an_evolution.jl")
include("co2_concentration.jl")
include("leaf_boundary_resistance.jl")
include("photosynthesis_biochemical.jl")
include("soil_resistance.jl")
include("enhancement_factor_ra_pleim.jl")

end
