using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: lwr_abs_indoors_no_int_mass

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "BuildingEnergyModel.LWRabsIndoorsNoIntMass.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

@testset "MATLAB" begin
    LWRinB, LWRoutB, LWRabsB, LWREBB = lwr_abs_indoors_no_int_mass(
        input_vars["Tinwallsun"],
        input_vars["Tinwallshd"],
        input_vars["Tceiling"],
        input_vars["Tground"],
        input_vars["Hbuild"],
        input_vars["Wroof"],
        input_vars["ec"],
        input_vars["eg"],
        input_vars["ew"],
    )

    # Test incoming radiation
    @test LWRinB.LWRinCeiling ≈ output_vars["LWRinB"]["LWRinCeiling"]
    @test LWRinB.LWRinWallsun ≈ output_vars["LWRinB"]["LWRinWallsun"]
    @test LWRinB.LWRinWallshd ≈ output_vars["LWRinB"]["LWRinWallshd"]
    @test LWRinB.LWRinGround ≈ output_vars["LWRinB"]["LWRinGround"]

    # Test outgoing radiation
    @test LWRoutB.LWRoutCeiling ≈ output_vars["LWRoutB"]["LWRoutCeiling"]
    @test LWRoutB.LWRoutWallsun ≈ output_vars["LWRoutB"]["LWRoutWallsun"]
    @test LWRoutB.LWRoutWallshd ≈ output_vars["LWRoutB"]["LWRoutWallshd"]
    @test LWRoutB.LWRoutGround ≈ output_vars["LWRoutB"]["LWRoutGround"]

    # Test absorbed radiation
    @test LWRabsB.LWRabsCeiling ≈ output_vars["LWRabsB"]["LWRabsCeiling"] atol=1e-11
    @test LWRabsB.LWRabsWallsun ≈ output_vars["LWRabsB"]["LWRabsWallsun"] atol=1e-11
    @test LWRabsB.LWRabsWallshd ≈ output_vars["LWRabsB"]["LWRabsWallshd"] atol=1e-11
    @test LWRabsB.LWRabsGround ≈ output_vars["LWRabsB"]["LWRabsGround"] atol=1e-11

    # Test energy balance
    @test LWREBB.LWREBCeiling ≈ output_vars["LWREBB"]["LWREBCeiling"] atol=1e-11
    @test LWREBB.LWREBWallsun ≈ output_vars["LWREBB"]["LWREBWallsun"] atol=1e-11
    @test LWREBB.LWREBWallshd ≈ output_vars["LWREBB"]["LWREBWallshd"] atol=1e-11
    @test LWREBB.LWREBGround ≈ output_vars["LWREBB"]["LWREBGround"] atol=1e-11
end
