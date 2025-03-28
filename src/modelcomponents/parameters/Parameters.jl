module Parameters

using TethysChlorisCore
using TethysChlorisCore.ModelComponents

include("SoilParameters.jl")
include("SurfaceFractions.jl")
include("ThermalProperties.jl")
include("UrbanGeometryParameters.jl")
include("VegetationParameters.jl")

# Depends on SurfaceFractions
include("OpticalProperties.jl")

include("ParameterSet.jl")

export UrbanGeometryParameters

end
