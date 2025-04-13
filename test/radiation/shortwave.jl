using Test
using UrbanTethysChloris.Radiation: direct_shortwave_trees

FT = Float64

@testset "direct_shortwave_trees" begin
    @testset "MATLAB - completely shaded" begin
        # Both trees are completely shaded
        d_tree = FT(0.136505190311419)
        h_can = FT(1.124567474048443)
        h_tree = FT(0.815051903114187)
        r_tree = FT(0.05)
        SWR_dir = FT(0.0)
        theta_n = FT(-2.820553888182081)
        theta_Z = FT(pi / 2.0)

        SWR_tree1, SWR_tree2 = direct_shortwave_trees(
            h_can, d_tree, h_tree, r_tree, theta_Z, theta_n, SWR_dir
        )

        @test SWR_tree1 == 0.0
        @test SWR_tree2 == 0.0
    end

    @testset "MATLAB - partially sunlit tree1" begin
        # tree 1 is partially sunlit
        # tree 2 is completely shaded
        d_tree = FT(0.136505190311419)
        h_can = FT(1.124567474048443)
        h_tree = FT(0.815051903114187)
        r_tree = FT(0.05)
        SWR_dir = FT(1.027251526757701)
        theta_n = FT(-0.593239662921201)
        theta_Z = FT(1.356513615208954)

        SWR_tree1, SWR_tree2 = direct_shortwave_trees(
            h_can, d_tree, h_tree, r_tree, theta_Z, theta_n, SWR_dir
        )

        @test SWR_tree1 ≈ 0.674302217356485
        @test SWR_tree2 == 0.0
    end

    @testset "MATLAB - completely / partially sunlit" begin
        # tree 1 is completely sunlit
        # tree 2 is partially sunlit
        d_tree = FT(0.136505190311419)
        h_can = FT(1.124567474048443)
        h_tree = FT(0.815051903114187)
        r_tree = FT(0.05)
        SWR_dir = FT(1.832071973971378)
        theta_n = FT(-0.109441709557626)
        theta_Z = FT(1.234325341923081)

        SWR_tree1, SWR_tree2 = direct_shortwave_trees(
            h_can, d_tree, h_tree, r_tree, theta_Z, theta_n, SWR_dir
        )

        @test SWR_tree1 ≈ 0.610938531288319
        @test SWR_tree2 ≈ 0.537875671264354
    end

    @testset "MATLAB - completely sunlit" begin
        # both trees are completely sunlit
        d_tree = FT(0.222413793103448)
        h_can = FT(1.120689655172414)
        h_tree = FT(0.812068965517241)
        r_tree = FT(0.05)
        SWR_dir = FT(350.4668408605486)
        theta_n = FT(2.953263261382311)
        theta_Z = FT(1.051771637797578)

        SWR_tree1, SWR_tree2 = direct_shortwave_trees(
            h_can, d_tree, h_tree, r_tree, theta_Z, theta_n, SWR_dir
        )

        @test SWR_tree1 ≈ 117.3950854817185
        @test SWR_tree2 ≈ 117.3950854817185
    end
end
