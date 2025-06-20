module ConductiveHeat

using ..ModelComponents
using ..Soil: soil_parameters_total, soil_thermal_properties
using LinearAlgebra: dot

include("conductive_heat_flux_green_roof.jl")
include("conductive_heat_flux_roof_imp.jl")

end
