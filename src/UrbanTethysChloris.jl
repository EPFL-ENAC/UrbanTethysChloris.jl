module UrbanTethysChloris

using TethysChlorisCore
using Dates
using Statistics

include(joinpath("modelcomponents", "ModelComponents.jl"))
using .ModelComponents

include(joinpath("ray_tracing", "RayTracing.jl"))
using .RayTracing

include(joinpath("radiation", "Radiation.jl"))
using .Radiation

include("incoming_longwave.jl")
include("set_sun_variables.jl")

end
