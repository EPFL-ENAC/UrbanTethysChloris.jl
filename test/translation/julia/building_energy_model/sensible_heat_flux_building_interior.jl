using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: sensible_heat_flux_building_interior

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "BuildingEnergyModel.SensibleHeatFluxBuildingInterior.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

@testset "MATLAB" begin
    HbinWallSun, HbinWallshd, HBinRoof, HBinGround, HbinIntMass, HbinWindow = sensible_heat_flux_building_interior(
        input_vars["Tbin"],
        input_vars["Tinwallsun"],
        input_vars["Tinwallshd"],
        input_vars["Tceiling"],
        input_vars["Tinground"],
        input_vars["Tintmass"],
        input_vars["Twindow"],
    )

    @test HbinWallSun ≈ output_vars["HbinWallSun"]
    @test HbinWallshd ≈ output_vars["HbinWallshd"]
    @test HBinRoof ≈ output_vars["HBinRoof"]
    @test HBinGround ≈ output_vars["HBinGround"]
    @test HbinIntMass ≈ output_vars["HbinIntMass"]
    @test HbinWindow ≈ output_vars["HbinWindow"]
end
