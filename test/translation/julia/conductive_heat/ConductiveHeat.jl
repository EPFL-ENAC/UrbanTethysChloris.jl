using SafeTestSets

@safetestset "ConductiveHeatFlux_GreenRoof" begin
    include("conductive_heat_flux_green_roof.jl")
end

@safetestset "ConductiveHeatFlux_RoofImp" begin
    include("conductive_heat_flux_roof_imp.jl")
end

@safetestset "ConductiveHeatFlux_Walls" begin
    include("conductive_heat_flux_walls.jl")
end

@safetestset "ConductiveHeatFluxFR_GroundImp" begin
    include("conductive_heat_flux_ground_fr.jl")
end

@safetestset "ConductiveHeatFluxFR_GroundVegBare" begin
    include("conductive_heat_flux_ground_vb.jl")
end

@safetestset "Soil_Heat" begin
    include("soil_heat.jl")
end
