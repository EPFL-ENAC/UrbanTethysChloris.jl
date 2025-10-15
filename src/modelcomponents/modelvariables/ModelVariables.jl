module ModelVariables

using TethysChlorisCore
using TethysChlorisCore.ModelComponents

abstract type AbstractModelVariables{FT<:AbstractFloat} <:
              AbstractIndividualModelComponent{FT} end

include("SolverVariables.jl")

end
