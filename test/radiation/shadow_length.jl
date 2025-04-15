using Test
using UrbanTethysChloris.Radiation: shadow_length_no_tree, shadow_length_with_trees

FT = Float64

@testset "shadow_length_no_tree" begin
    @testset "MATLAB" begin
        # Test case parameters
        h_can = FT(1.124567474048443)
        theta_n = FT(-2.820553888182081)
        theta_Z = FT(pi / 2.0)
        w_can = FT(1.0)

        X_shadow, X_tree, n_shadow, n_tree = shadow_length_no_tree(
            h_can, w_can, theta_Z, theta_n
        )

        @test X_shadow ≈ 1.0
        @test X_tree ≈ 0.0
        @test n_shadow ≈ 1.0
        @test n_tree ≈ 0.0
    end
end

@testset "shadow_length_with_trees" begin
    @testset "MATLAB" begin
        # Test case parameters
        h_can = FT(1.124567474048443)
        w_can = FT(1.0)
        d_tree = FT(0.136505190311419)
        h_tree = FT(0.815051903114187)
        r_tree = FT(0.05)
        theta_n = FT(-2.820553888182081)
        theta_Z = FT(pi / 2.0)

        X_shadow, X_tree, n_shadow, n_tree = shadow_length_with_trees(
            h_can, w_can, d_tree, h_tree, r_tree, theta_Z, theta_n
        )

        @test X_shadow ≈ 1.0
        @test X_tree ≈ 0.0
        @test n_shadow ≈ 1.0
        @test n_tree ≈ 0.0
    end
end
