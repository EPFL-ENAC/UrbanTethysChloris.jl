module Parameters

using TethysChlorisCore
using TethysChlorisCore.ModelComponents

include("PersonParameters.jl")
include("SurfaceFractions.jl")
export LocationSpecificSurfaceFractions
include("ThermalProperties.jl")
export LocationSpecificThermalProperties
include("UrbanGeometryParameters.jl")
export UrbanGeometryParameters
include("VegetationParameters.jl")
export HeightDependentVegetationParameters
include("SoilParameters.jl")
export VegetatedSoilParameters
include("BuildingEnergyModelParameters.jl")
export IndoorOpticalProperties, HVACParameters, ThermalBuilding, WindowParameters
include("LocationProperties.jl")

# Depends on SurfaceFractions
include("OpticalProperties.jl")
export VegetatedOpticalProperties, SimpleOpticalProperties

include("ParameterSet.jl")

export initialize_parameter_set

end
