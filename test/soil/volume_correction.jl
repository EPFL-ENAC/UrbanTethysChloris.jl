using Test
using UrbanTethysChloris.Soil: volume_correction

@testset "MATLAB" begin
    @testset "Base case" begin
        # Test setup
        V = [-0.1, 0.5, -0.2, 0.3, -0.15]
        EvL_Zs = [-1.0, 1.0, 0.0, 0.0, 0.0]
        RfH_Zs = [1.0 0.5; 1.0 1.0; 1.0 0.5; 0.0 0.0; 0.0 0.0]
        RfL_Zs = [-1.0 -0.5; 1.0 0.8; 0.0 0.2; 0.0 0.0; 0.0 0.0]
        EG = 0.2
        T_H = [0.3, 0.2, 0.1]
        T_L = [0.25, 0.15, 0.2]
        Lk = 0.1

        # Run correction
        V_new, T_H_new, T_L_new, EG_new, Lk_new = volume_correction(
            V, EvL_Zs, RfH_Zs, RfL_Zs, EG, T_H, T_L, Lk
        )

        # Test results
        @test all(V_new .== [0.0, 0.2, 0.0, 0.15, 0.0])
        @test all(T_H_new .== [0.3, 0.2, 0.1])
        @test all(T_L_new .== [0.25, 0.15, 0.2])
        @test Lk_new == Lk
        @test EG_new == EG
    end

    @testset "L correction" begin
        EG = 0.0
        EvL_Zs = [1.0, 0.0, 0.0, 0.0]
        Lk = 0.0
        RfH_Zs = reshape([0.0, 0.0, 0.0, 0.0], 1, 4)
        RfL_Zs = reshape(
            [0.284979001762539, 0.207810379282956, 0.342622081539846, 0.164602146793880],
            1,
            4,
        )
        T_H = [0.0, 0.0, 0.0, 0.0]
        T_L = [0.0, 0.0, 0.0, 0.0]
        V = [-0.961232650478129, -0.961232642886289, -2.883697928658866, -4.806163206839603]

        V_new, T_H_new, T_L_new, EG_new, Lk_new = volume_correction(
            V, EvL_Zs, RfH_Zs, RfL_Zs, EG, T_H, T_L, Lk
        )

        @test all(V_new .== 0)
        @test all(T_H_new .== 0)
        @test all(T_L_new .== 0)
        @test EG_new ≈ -9.612326428862886
        @test Lk_new == 0
    end
end
