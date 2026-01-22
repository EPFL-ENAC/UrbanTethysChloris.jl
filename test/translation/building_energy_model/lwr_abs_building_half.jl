using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: lwr_abs_building_half
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("BuildingEnergyModel.LWRabsBuildingHalf.json")

@testset "MATLAB" begin
    LWRinB, LWRoutB, LWRabsB = lwr_abs_building_half(
        input_vars["Tceiling"],
        input_vars["Tinwall"],
        input_vars["Tintmass"],
        input_vars["Tground"],
        input_vars["A_c"],
        input_vars["A_g"],
        input_vars["A_h"],
        input_vars["F_gc"],
        input_vars["F_gw"],
        input_vars["F_ww"],
        input_vars["F_wg"],
        input_vars["F_wc"],
        input_vars["F_cg"],
        input_vars["F_cw"],
        input_vars["ec"],
        input_vars["eg"],
        input_vars["ew"],
        input_vars["em"],
    )

    @test LWRinB.LWRinCeiling ≈ output_vars["LWRinB"]["LWRinCeiling"]
    @test LWRinB.LWRinWall ≈ output_vars["LWRinB"]["LWRinWall"]
    @test LWRinB.LWRinInternalMass ≈ output_vars["LWRinB"]["LWRinInternalMass"]
    @test LWRinB.LWRinGround ≈ output_vars["LWRinB"]["LWRinGround"]

    @test LWRoutB.LWRoutCeiling ≈ output_vars["LWRoutB"]["LWRoutCeiling"]
    @test LWRoutB.LWRoutWall ≈ output_vars["LWRoutB"]["LWRoutWall"]
    @test LWRoutB.LWRoutInternalMass ≈ output_vars["LWRoutB"]["LWRoutInternalMass"]
    @test LWRoutB.LWRoutGround ≈ output_vars["LWRoutB"]["LWRoutGround"]

    @test LWRabsB.LWRabsCeiling ≈ output_vars["LWRabsB"]["LWRabsCeiling"]
    @test LWRabsB.LWRabsWall ≈ output_vars["LWRabsB"]["LWRabsWall"]
    @test LWRabsB.LWRabsInternalMass ≈ output_vars["LWRabsB"]["LWRabsInternalMass"]
    @test LWRabsB.LWRabsGround ≈ output_vars["LWRabsB"]["LWRabsGround"]
end
