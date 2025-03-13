module Parameters

using TethysChlorisCore
using TethysChlorisCore.ModelComponents

include("UrbanGeometryParameters.jl")
include("SurfaceFractions.jl")
include("ThermalProperties.jl")
include("VegetationParameters.jl")

include("ParameterSet.jl")

export UrbanGeometryParameters

end
