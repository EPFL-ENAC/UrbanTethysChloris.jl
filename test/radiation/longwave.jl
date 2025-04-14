using Test
using UrbanTethysChloris.Radiation: longwave_absorbed_no_tree, longwave_absorbed_with_trees
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

@testset "longwave_absorbed_with_trees" begin
    @testset "MATLAB" begin
        egbare = FT(0.95)
        egimp = FT(0.95)
        egveg = FT(0.917997045345903)
        et = FT(0.994483435579239)
        ew = FT(0.95)
        fgbare = FT(0.0)
        fgimp = FT(0.455)
        fgveg = FT(0.545)
        h_can = FT(1.124567474048443)
        LWR = FT(320.4495620328094)
        r_tree = FT(0.05)
        Tgbare = Tgimp = Tgveg = Ttree = Twshade = Twsun = FT(277.05)
        w_can = FT(1.0)
        view_factor = ViewFactor{FT}(
            F_gs_T=0.320368159203980,
            F_gt_T=0.097393034825871,
            F_gw_T=0.291119402985075,
            F_ww_T=0.361365059318791,
            F_wt_T=0.150343283582090,
            F_wg_T=0.258872330654420,
            F_ws_T=0.229419326444700,
            F_sg_T=0.320368159203980,
            F_sw_T=0.257997512437811,
            F_st_T=0.163636815920398,
            F_tg_T=0.155005829152584,
            F_tw_T=0.269085119079438,
            F_ts_T=0.260436081255499,
            F_tt_T=0.046387851433041,
        )

        LWRin_nT, LWRout_nT, LWRabs_nT, LWREB_nT = longwave_absorbed_with_trees(
            h_can,
            w_can,
            r_tree,
            LWR,
            fgveg,
            fgbare,
            fgimp,
            ew,
            et,
            egveg,
            egbare,
            egimp,
            Tgimp,
            Tgbare,
            Tgveg,
            Twsun,
            Twshade,
            Ttree,
            view_factor,
        )

        @test LWRin_nT.GroundImp ≈ 329.598147532779
        @test LWRin_nT.GroundBare == 0
        @test LWRin_nT.GroundVeg ≈ 329.598147532779
        @test LWRin_nT.Tree ≈ 330.3750905769887
        @test LWRin_nT.WallSun ≈ 330.7925407436207
        @test LWRin_nT.WallShade ≈ 330.7925407436207
        @test LWRin_nT.TotalGround ≈ 329.598147532779
        @test LWRin_nT.TotalCanyon ≈ 1281.176002986166

        @test LWRout_nT.GroundImp ≈ 333.8304769888871
        @test LWRout_nT.GroundBare == 0
        @test LWRout_nT.GroundVeg ≈ 333.6879011492398
        @test LWRout_nT.Tree ≈ 334.0329404712533
        @test LWRout_nT.WallSun ≈ 333.8901966494292
        @test LWRout_nT.WallShade ≈ 333.8901966494292
        @test LWRout_nT.TotalGround ≈ 333.7527731562794
        @test LWRout_nT.TotalCanyon ≈ 1294.595969635749

        @test LWRabs_nT.GroundImp ≈ -4.232329456108342
        @test LWRabs_nT.GroundBare == 0
        @test LWRabs_nT.GroundVeg ≈ -4.089753616459356
        @test LWRabs_nT.Tree ≈ -3.657849894255961
        @test LWRabs_nT.WallSun ≈ -3.097655905808094
        @test LWRabs_nT.WallShade ≈ -3.097655905806957
        @test LWRabs_nT.TotalGround ≈ -4.154625623499645
        @test LWRabs_nT.TotalCanyon ≈ -13.419966649575812
    end
end
