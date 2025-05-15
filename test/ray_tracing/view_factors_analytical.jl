using Test
using UrbanTethysChloris.RayTracing: view_factors_analytical

FT = Float64

@testset "view_factors_analytical" begin
    @testset "MATLAB" begin
        H = FT(6.5)
        W = FT(5.8)

        vf = view_factors_analytical(H, W)

        @test vf.F_gs ≈ 0.381290805156702
        @test vf.F_gw ≈ 0.309354597421649
        @test vf.F_ww ≈ 0.447921026139827
        @test vf.F_wg ≈ 0.276039486930087
        @test vf.F_ws ≈ vf.F_wg
        @test vf.F_sg ≈ vf.F_gs
        @test vf.F_sw ≈ vf.F_gw
    end
end
