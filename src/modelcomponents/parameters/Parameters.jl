module Parameters

using TethysChlorisCore
using TethysChlorisCore.ModelComponents
import ...UrbanTethysChloris: check_extraneous_fields

include("PersonParameters.jl")
include("SoilParameters.jl")
include("SurfaceFractions.jl")
include("ThermalProperties.jl")
include("UrbanGeometryParameters.jl")
include("VegetationParameters.jl")
include("BuildingEnergyModelParameters.jl")

# Depends on SurfaceFractions
include("OpticalProperties.jl")

include("ParameterSet.jl")

export UrbanGeometryParameters

end
