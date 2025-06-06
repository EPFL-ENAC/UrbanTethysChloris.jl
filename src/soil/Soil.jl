module Soil

include("conductivity_suction.jl")
include("evaporation_layers.jl")
include("infiltration.jl")
include("leakage_bottom.jl")
include("root_fraction_general.jl")
include("root_soil_conductance.jl")
include("soil_moistures_rich_comp.jl")
include("soil_moistures_rich_comp_lat2.jl")
include("soil_moistures_rich_comp_lat3.jl")
include("soil_parameters.jl")
include("soil_parameters2.jl")
include("soil_parameters_total.jl")
include("soil_thermal_properties.jl")
include("soil_water_multilayer.jl")
include("volume_correction.jl")
include("soil_moisture_conductivity_update.jl")

end
