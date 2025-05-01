module ModelComponents

include(joinpath("parameters", "Parameters.jl"))
using .Parameters

include(joinpath("forcinginputs", "ForcingInputs.jl"))
using .ForcingInputs

end
