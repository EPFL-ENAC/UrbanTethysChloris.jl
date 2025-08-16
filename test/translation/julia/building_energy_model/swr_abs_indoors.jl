using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: swr_abs_indoors
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("BuildingEnergyModel.SWRabsIndoors.mat")

@testset "MATLAB" begin
    SWRinB, SWRoutB, SWRabsB = swr_abs_indoors(
        input_vars["SWRinWsun"],
        input_vars["SWRinWshd"],
        input_vars["Hbuild"],
        input_vars["Wroof"],
        input_vars["abc"],
        input_vars["abw"],
        input_vars["abg"],
        input_vars["abm"],
    )

    @test SWRinB.SWRinCeiling ≈ output_vars["SWRinB"]["SWRinCeiling"]
    @test SWRinB.SWRinWallsun ≈ output_vars["SWRinB"]["SWRinWallsun"]
    @test SWRinB.SWRinWallshd ≈ output_vars["SWRinB"]["SWRinWallshd"]
    @test SWRinB.SWRinInternalMass ≈ output_vars["SWRinB"]["SWRinInternalMass"]
    @test SWRinB.SWRinGround ≈ output_vars["SWRinB"]["SWRinGround"]

    @test SWRoutB.SWRoutCeiling ≈ output_vars["SWRoutB"]["SWRoutCeiling"]
    @test SWRoutB.SWRoutWallsun ≈ output_vars["SWRoutB"]["SWRoutWallsun"]
    @test SWRoutB.SWRoutWallshd ≈ output_vars["SWRoutB"]["SWRoutWallshd"]
    @test SWRoutB.SWRoutInternalMass ≈ output_vars["SWRoutB"]["SWRoutInternalMass"]
    @test SWRoutB.SWRoutGround ≈ output_vars["SWRoutB"]["SWRoutGround"]

    @test SWRabsB.SWRabsCeiling ≈ output_vars["SWRabsB"]["SWRabsCeiling"]
    @test SWRabsB.SWRabsWallsun ≈ output_vars["SWRabsB"]["SWRabsWallsun"]
    @test SWRabsB.SWRabsWallshd ≈ output_vars["SWRabsB"]["SWRabsWallshd"]
    @test SWRabsB.SWRabsInternalMass ≈ output_vars["SWRabsB"]["SWRabsInternalMass"]
    @test SWRabsB.SWRabsGround ≈ output_vars["SWRabsB"]["SWRabsGround"]
end
