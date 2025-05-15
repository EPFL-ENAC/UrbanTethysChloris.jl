module Parameters

using TethysChlorisCore
using TethysChlorisCore.ModelComponents

include("PersonParameters.jl")
include("SoilParameters.jl")
include("SurfaceFractions.jl")
include("ThermalProperties.jl")
include("UrbanGeometryParameters.jl")
include("VegetationParameters.jl")
include("BuildingEnergyModelParameters.jl")
include("LocationProperties.jl")

# Depends on SurfaceFractions
include("OpticalProperties.jl")

include("ParameterSet.jl")

export UrbanGeometryParameters

end
