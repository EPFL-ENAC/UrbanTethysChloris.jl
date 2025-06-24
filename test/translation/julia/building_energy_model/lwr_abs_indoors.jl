using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: lwr_abs_indoors

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "building_energy_model.LWRabsIndoors.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

@testset "MATLAB" begin
    LWRinB, LWRoutB, LWRabsB = lwr_abs_indoors(
        input_vars["Tinwallsun"],
        input_vars["Tinwallshd"],
        input_vars["Tceiling"],
        input_vars["Tground"],
        input_vars["Tintmass"],
        input_vars["Hbuild"],
        input_vars["Wroof"],
        input_vars["ec"],
        input_vars["eg"],
        input_vars["em"],
        input_vars["ew"],
    )

    # Test incoming radiation
    @test LWRinB.LWRinCeiling ≈ output_vars["LWRinB"]["LWRinCeiling"]
    @test LWRinB.LWRinWallsun ≈ output_vars["LWRinB"]["LWRinWallsun"]
    @test LWRinB.LWRinWallshd ≈ output_vars["LWRinB"]["LWRinWallshd"]
    @test LWRinB.LWRinInternalMass ≈ output_vars["LWRinB"]["LWRinInternalMass"]
    @test LWRinB.LWRinGround ≈ output_vars["LWRinB"]["LWRinGround"]

    # Test outgoing radiation
    @test LWRoutB.LWRoutCeiling ≈ output_vars["LWRoutB"]["LWRoutCeiling"]
    @test LWRoutB.LWRoutWallsun ≈ output_vars["LWRoutB"]["LWRoutWallsun"]
    @test LWRoutB.LWRoutWallshd ≈ output_vars["LWRoutB"]["LWRoutWallshd"]
    @test LWRoutB.LWRoutInternalMass ≈ output_vars["LWRoutB"]["LWRoutInternalMass"]
    @test LWRoutB.LWRoutGround ≈ output_vars["LWRoutB"]["LWRoutGround"]

    # Test absorbed radiation
    @test LWRabsB.LWRabsCeiling ≈ output_vars["LWRabsB"]["LWRabsCeiling"]
    @test LWRabsB.LWRabsWallsun ≈ output_vars["LWRabsB"]["LWRabsWallsun"]
    @test LWRabsB.LWRabsWallshd ≈ output_vars["LWRabsB"]["LWRabsWallshd"]
    @test LWRabsB.LWRabsInternalMass ≈ output_vars["LWRabsB"]["LWRabsInternalMass"]
    @test LWRabsB.LWRabsGround ≈ output_vars["LWRabsB"]["LWRabsGround"]
end
