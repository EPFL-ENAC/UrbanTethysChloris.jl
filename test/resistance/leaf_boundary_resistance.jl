using Test
using UrbanTethysChloris.Resistance: leaf_boundary_resistance

FT = Float64

@testset "MATLAB" begin
    Ta = 4.049999999999955
    Ts = 3.899999999999977

    @testset "High-level" begin
        d_leaf = 0.8
        disp_h = 0.1
        hc = 0.15
        LAI = 2.5

        Ws = 0.057707208970656
        zatm = 25
        zom = 0.01845

        rb = leaf_boundary_resistance(Ws, Ts, Ta, hc, d_leaf, LAI, zatm, disp_h, zom)
        @test rb ≈ 200.5974949314946
    end

    @testset "Low-level" begin
        alpha = 0.350645182969221
        d_leaf = 4.0
        u_hc = 0.133882970988947

        rb = leaf_boundary_resistance(u_hc, Ts, Ta, d_leaf, alpha)
        @test rb ≈ 119.182380966834
    end
end
