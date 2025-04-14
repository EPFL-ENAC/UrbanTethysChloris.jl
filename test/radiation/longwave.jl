using Test
using UrbanTethysChloris.Radiation: longwave_absorbed_no_tree
using UrbanTethysChloris.RayTracing: ViewFactor

FT = Float64

@testset "longwave_absorbed_no_tree" begin
    @testset "MATLAB" begin
        egbare = FT(0.95)
        egimp = FT(0.95)
        egveg = FT(0.917997045345903)
        ew = FT(0.95)
        fgbare = FT(0.0)
        fgimp = FT(0.455)
        fgveg = FT(0.545)
        h_can = FT(1.124567474048443)
        LWR = FT(320.4495620328094)
        Tgbare = Tgimp = Tgveg = Twshade = Twsun = FT(277.05)
        w_can = FT(1.0)
        view_factor = ViewFactor{FT}(
            F_gs_nT=0.380308601808867,
            F_gw_nT=0.309845699095567,
            F_ww_nT=0.448951341300808,
            F_wg_nT=0.275524329349596,
            F_ws_nT=0.275524329349596,
            F_sg_nT=0.380308601808867,
            F_sw_nT=0.309845699095567,
        )

        LWRin_nT, LWRout_nT, LWRabs_nT, LWREB_nT = longwave_absorbed_no_tree(
            h_can,
            w_can,
            LWR,
            fgveg,
            fgbare,
            fgimp,
            ew,
            egveg,
            egbare,
            egimp,
            Tgimp,
            Tgbare,
            Tgveg,
            Twsun,
            Twshade,
            view_factor,
        )

        @test LWRin_nT.GroundImp ≈ 328.7577185187801
        @test LWRin_nT.GroundBare == 0
        @test LWRin_nT.GroundVeg ≈ 328.7577185187801
        @test LWRin_nT.Tree == 0
        @test LWRin_nT.WallSun ≈ 330.1183609102541
        @test LWRin_nT.WallShade ≈ 330.1183609102541
        @test LWRin_nT.TotalGround ≈ 328.7577185187801
        @test LWRin_nT.TotalCanyon ≈ 1071.238461050493

        @test LWRout_nT.GroundImp ≈ 333.7884555381871
        @test LWRout_nT.GroundBare == 0
        @test LWRout_nT.GroundVeg ≈ 333.6189834869148
        @test LWRout_nT.Tree == 0
        @test LWRout_nT.WallSun ≈ 333.8564876577608
        @test LWRout_nT.WallShade ≈ 333.8564876577608
        @test LWRout_nT.TotalGround ≈ 333.6960932702438
        @test LWRout_nT.TotalCanyon ≈ 1084.58438731019

        @test LWRabs_nT.GroundImp ≈ -5.030737019408202
        @test LWRabs_nT.GroundBare == 0
        @test LWRabs_nT.GroundVeg ≈ -4.861264968134406
        @test LWRabs_nT.Tree == 0
        @test LWRabs_nT.WallSun ≈ -3.738126747507519
        @test LWRabs_nT.WallShade ≈ -3.738126747507519
        @test LWRabs_nT.TotalGround ≈ -4.938374751463984
        @test LWRabs_nT.TotalCanyon ≈ -13.345926259698887
    end
end
