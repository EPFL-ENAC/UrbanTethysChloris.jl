using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: swr_abs_indoors_no_int_mass

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "BuildingEnergyModel.SWRabsIndoorsNoIntMass.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

@testset "MATLAB" begin
    SWRinB, SWRoutB, SWRabsB, SWREBB = swr_abs_indoors_no_int_mass(
        input_vars["SWRinWsun"],
        input_vars["SWRinWshd"],
        input_vars["Hbuild"],
        input_vars["Wroof"],
        input_vars["abc"],
        input_vars["abw"],
        input_vars["abg"],
    )

    # Test incoming radiation
    @test SWRinB.SWRinCeiling ≈ output_vars["SWRinB"]["SWRinCeiling"]
    @test SWRinB.SWRinWallsun ≈ output_vars["SWRinB"]["SWRinWallsun"]
    @test SWRinB.SWRinWallshd ≈ output_vars["SWRinB"]["SWRinWallshd"]
    @test SWRinB.SWRinGround ≈ output_vars["SWRinB"]["SWRinGround"]

    # Test outgoing radiation
    @test SWRoutB.SWRoutCeiling ≈ output_vars["SWRoutB"]["SWRoutCeiling"]
    @test SWRoutB.SWRoutWallsun ≈ output_vars["SWRoutB"]["SWRoutWallsun"]
    @test SWRoutB.SWRoutWallshd ≈ output_vars["SWRoutB"]["SWRoutWallshd"]
    @test SWRoutB.SWRoutGround ≈ output_vars["SWRoutB"]["SWRoutGround"]

    # Test absorbed radiation
    @test SWRabsB.SWRabsCeiling ≈ output_vars["SWRabsB"]["SWRabsCeiling"]
    @test SWRabsB.SWRabsWallsun ≈ output_vars["SWRabsB"]["SWRabsWallsun"]
    @test SWRabsB.SWRabsWallshd ≈ output_vars["SWRabsB"]["SWRabsWallshd"]
    @test SWRabsB.SWRabsGround ≈ output_vars["SWRabsB"]["SWRabsGround"]

    # Test energy balance
    @test SWREBB.SWREBCeiling ≈ output_vars["SWREBB"]["SWREBCeiling"]
    @test SWREBB.SWREBWallsun ≈ output_vars["SWREBB"]["SWREBWallsun"]
    @test SWREBB.SWREBWallshd ≈ output_vars["SWREBB"]["SWREBWallshd"]
    @test SWREBB.SWREBGround ≈ output_vars["SWREBB"]["SWREBGround"]
end
