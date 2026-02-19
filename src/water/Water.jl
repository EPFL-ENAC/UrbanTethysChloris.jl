module Water

using OrdinaryDiffEqRosenbrock: ODEProblem, Rosenbrock23, solve
using ADTypes: AutoFiniteDiff
using TethysChlorisCore: AbstractODEOptions, ODEOptions
using ..ModelComponents
using ..Soil

include("water_canyon.jl")
include("water_ground.jl")
include("water_impervious.jl")
include("water_roof.jl")
include("water_soil.jl")
include("water_vegetation.jl")

end
