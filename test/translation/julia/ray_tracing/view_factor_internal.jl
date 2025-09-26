using Test
using MAT
using UrbanTethysChloris.RayTracing: view_factor_internal
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("BuildingEnergyModel.ViewFactorInternal.json")

@testset "MATLAB" begin
    F_gc, F_gw, F_ww, F_wg, F_wc, F_cg, F_cw, ViewFactor = view_factor_internal(
        input_vars["Hbuild"], input_vars["Wroof"]
    )

    @test ViewFactor.F_gc ≈ output_vars["ViewFactor"]["F_gc"]
    @test ViewFactor.F_gw ≈ output_vars["ViewFactor"]["F_gw"]
    @test ViewFactor.F_ww ≈ output_vars["ViewFactor"]["F_ww"]
    @test ViewFactor.F_wg ≈ output_vars["ViewFactor"]["F_wg"]
    @test ViewFactor.F_wc ≈ output_vars["ViewFactor"]["F_wc"]
    @test ViewFactor.F_cg ≈ output_vars["ViewFactor"]["F_cg"]
    @test ViewFactor.F_cw ≈ output_vars["ViewFactor"]["F_cw"]
end
