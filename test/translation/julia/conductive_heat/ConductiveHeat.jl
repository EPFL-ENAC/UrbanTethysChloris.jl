using SafeTestSets

@safetestset "ConductiveHeatFlux_GreenRoof" begin
    include("conductive_heat_flux_green_roof.jl")
end

@safetestset "ConductiveHeatFlux_RoofImp" begin
    include("conductive_heat_flux_roof_imp.jl")
end
