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

@safetestset "LWRabsBuildingHalf" begin
    include("lwr_abs_building_half.jl")
end

@safetestset "LWRabsIndoors" begin
    include("lwr_abs_indoors.jl")
end

@safetestset "LWRabsIndoorsNoIntMass" begin
    include("lwr_abs_indoors_no_int_mass.jl")
end

@safetestset "SWRabsIndoorsNoIntMass" begin
    include("swr_abs_indoors_no_int_mass.jl")
end

@safetestset "SensibleHeatFluxBuildingInterior" begin
    include("sensible_heat_flux_building_interior.jl")
end
