using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: swr_abs_building_half
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("BuildingEnergyModel.SWRabsBuildingHalf.mat")

@testset "MATLAB" begin
    SWRinB, SWRoutB, SWRabsB = swr_abs_building_half(
        input_vars["A_c"],
        input_vars["A_g"],
        input_vars["A_h"],
        input_vars["SWRinW"],
        input_vars["F_gc"],
        input_vars["F_gw"],
        input_vars["F_ww"],
        input_vars["F_wg"],
        input_vars["F_wc"],
        input_vars["F_cg"],
        input_vars["F_cw"],
        input_vars["abc"],
        input_vars["abw"],
        input_vars["abm"],
        input_vars["abg"],
    )

    @test SWRinB.SWRinCeiling ≈ output_vars["SWRinB"]["SWRinCeiling"]
    @test SWRinB.SWRinWall ≈ output_vars["SWRinB"]["SWRinWall"]
    @test SWRinB.SWRinInternalMass ≈ output_vars["SWRinB"]["SWRinInternalMass"]
    @test SWRinB.SWRinGround ≈ output_vars["SWRinB"]["SWRinGround"]

    @test SWRoutB.SWRoutCeiling ≈ output_vars["SWRoutB"]["SWRoutCeiling"]
    @test SWRoutB.SWRoutWall ≈ output_vars["SWRoutB"]["SWRoutWall"]
    @test SWRoutB.SWRoutInternalMass ≈ output_vars["SWRoutB"]["SWRoutInternalMass"]
    @test SWRoutB.SWRoutGround ≈ output_vars["SWRoutB"]["SWRoutGround"]

    @test SWRabsB.SWRabsCeiling ≈ output_vars["SWRabsB"]["SWRabsCeiling"]
    @test SWRabsB.SWRabsWall ≈ output_vars["SWRabsB"]["SWRabsWall"]
    @test SWRabsB.SWRabsInternalMass ≈ output_vars["SWRabsB"]["SWRabsInternalMass"]
    @test SWRabsB.SWRabsGround ≈ output_vars["SWRabsB"]["SWRabsGround"]
end
