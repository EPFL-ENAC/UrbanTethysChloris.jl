module Resistance

using ..Soil

include("aerodynamic_resistance.jl")
include("leaf_boundary_resistance.jl")
include("photosynthesis_biochemical.jl")
include("soil_resistance.jl")

end
