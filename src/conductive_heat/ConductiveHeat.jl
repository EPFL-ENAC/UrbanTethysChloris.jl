module ConductiveHeat

using ..ModelComponents
using ..Soil: soil_parameters_total, soil_thermal_properties
using LinearAlgebra: dot

export conductive_heat_flux_ground_vb

include("conductive_heat_flux_green_roof.jl")
include("conductive_heat_flux_roof_imp.jl")
include("conductive_heat_flux_walls.jl")
include("conductive_heat_flux_ground_imp.jl")
include("conductive_heat_flux_ground_fr.jl")
include("conductive_heat_flux_ground_vb.jl")
include("soil_heat.jl")

end
