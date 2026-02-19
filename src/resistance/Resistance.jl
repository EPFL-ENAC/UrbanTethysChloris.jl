module Resistance

using Roots
using ..ModelComponents
using ..RayTracing
using ..Soil
using ..Radiation
using ..UrbanTethysChloris: Model, ModelIttm
using SimpleNonlinearSolve: IntervalNonlinearProblem, solve, Brent
using SciMLBase: successful_retcode
using TethysChlorisCore: AbstractZeroFindingStrategies, AbstractModelOptions
using TethysChlorisCore: SimpleBrentStrategy, find_root

"""
    StomatalResistancePreCalc{FT<:AbstractFloat}

Precalculated stomatal resistance and CO2 concentration values for faster numerical solution.

This type stores precalculated stomatal resistance and internal CO2 concentration values
for both sunlit and shaded leaves, used in the energy balance and water flux calculations
for vegetation (roof, ground, and tree).

# Fields
- `rs_sun::FT`: Stomatal resistance for sunlit leaves [s/m]
- `rs_shd::FT`: Stomatal resistance for shaded leaves [s/m]
- `Ci_sun::FT`: Internal CO2 concentration for sunlit leaves [ppm]
- `Ci_shd::FT`: Internal CO2 concentration for shaded leaves [ppm]
"""
Base.@kwdef struct StomatalResistancePreCalc{FT<:AbstractFloat}
    rs_sun::FT
    rs_shd::FT
    Ci_sun::FT
    Ci_shd::FT
end

include("aerodynamic_resistance.jl")
include("backcalculate_obhukov_length.jl")
include("canopy_resistance_an_evolution.jl")
include("co2_concentration.jl")
include("enhancement_factor_ra_pleim.jl")
include("in_canyon_aerodynamic_resistance.jl")
include("leaf_boundary_resistance.jl")
include("photosynthesis_biochemical.jl")
include("precalculate_for_faster_numerical_solution.jl")
include("precalculate_stomatal_resistance_ground_tree.jl")
include("precalculate_stomatal_resistance_roof.jl")
include("soil_resistance.jl")
include("urban_roughness.jl")
include("wind_profile_canyon.jl")
include("wind_profile_point_output.jl")
include("wind_profile_roof.jl")

export StomatalResistancePreCalc

end
