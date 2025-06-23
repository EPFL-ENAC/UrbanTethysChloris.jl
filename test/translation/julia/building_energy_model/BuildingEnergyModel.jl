using SafeTestsets

@safetestset "ACHeatingModule" begin
    include("ac_heating_module.jl")
end

@safetestset "ACHeatingTurnOnOff" begin
    include("ac_heating_turn_on_off.jl")
end

@safetestset "ConductiveHeatFluxBuildingFloor" begin
    include("conductive_heat_flux_building_floor.jl")
end

# TODO: manually create test data
# @safetestset "LWRabsBuildingHalf" begin
#     include("lwr_abs_building_half.jl")
# end

# @safetestset "LWRabsIndoors" begin
#     include("lwr_abs_indoors.jl")
# end

@safetestset "LWRabsIndoorsNoIntMass" begin
    include("lwr_abs_indoors_no_int_mass.jl")
end

# @safetestset "SWRabsIndoors" begin
#     include("swr_abs_indoors.jl")
# end

@safetestset "SWRabsIndoorsNoIntMass" begin
    include("swr_abs_indoors_no_int_mass.jl")
end

@safetestset "SensibleHeatFluxBuildingInterior" begin
    include("sensible_heat_flux_building_interior.jl")
end

# @safetestset "SWRabsBuildingHalf" begin
#     include("swr_abs_building_half.jl")
# end

# @safetestset "HeatStorageChangeInternalMass" begin
#     include("heat_storage_change_internal_mass.jl")
# end

@safetestset "EBSolver_Building" begin
    include("eb_solver_building.jl")
end
