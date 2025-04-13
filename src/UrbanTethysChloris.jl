module UrbanTethysChloris

using TethysChlorisCore

include(joinpath("modelcomponents", "ModelComponents.jl"))
using .ModelComponents

include(joinpath("radiation", "Radiation.jl"))
using .Radiation

end
