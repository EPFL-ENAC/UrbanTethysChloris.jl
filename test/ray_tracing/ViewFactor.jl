using Test
using UrbanTethysChloris.RayTracing:
    ViewFactor, ViewFactorNoTrees, ViewFactorWithTrees, ViewFactorPoint

@testset "ViewFactor" begin
    @testset "Default constructors" begin
        vf_nt = ViewFactorNoTrees{Float64}()
        @test vf_nt.F_gs == 0.0
        @test vf_nt.F_gw == 0.0
        @test vf_nt.F_ww == 0.0

        vf_wt = ViewFactorWithTrees{Float64}()
        @test vf_wt.F_gs == 0.0
        @test vf_wt.F_gt == 0.0
        @test vf_wt.F_tt == 0.0
    end

    @testset "Combined constructor" begin
        vf_nt = ViewFactorNoTrees{Float64}(
            F_gs=0.5, F_gw=0.3, F_ww=0.1, F_wg=0.3, F_ws=0.6, F_sg=0.5, F_sw=0.5
        )

        vf_wt = ViewFactorWithTrees{Float64}(
            F_gs=0.4,
            F_gt=0.2,
            F_gw=0.2,
            F_ww=0.1,
            F_wt=0.2,
            F_wg=0.2,
            F_ws=0.5,
            F_sg=0.4,
            F_sw=0.3,
            F_st=0.3,
            F_tg=0.2,
            F_tw=0.2,
            F_ts=0.3,
            F_tt=0.3,
        )

        vf = ViewFactor(vf_nt, vf_wt)

        # Test transfer from no trees case
        @test vf.F_gs_nT == vf_nt.F_gs
        @test vf.F_gw_nT == vf_nt.F_gw
        @test vf.F_ww_nT == vf_nt.F_ww

        # Test transfer from with trees case
        @test vf.F_gs_T == vf_wt.F_gs
        @test vf.F_gt_T == vf_wt.F_gt
        @test vf.F_tt_T == vf_wt.F_tt
    end

    @testset "ViewFactorPoint" begin
        vfp = ViewFactorPoint{Float64}(0.2, 0.3, 0.1, 0.2, 0.2)
        @test vfp.F_pg == 0.2
        @test vfp.F_ps == 0.3
        @test vfp.F_pt == 0.1
        @test vfp.F_pwLeft == 0.2
        @test vfp.F_pwRight == 0.2
        @test sum([vfp.F_pg, vfp.F_ps, vfp.F_pt, vfp.F_pwLeft, vfp.F_pwRight]) ≈ 1.0
    end
end
