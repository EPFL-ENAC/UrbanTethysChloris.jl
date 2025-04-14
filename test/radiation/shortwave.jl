using Test
using UrbanTethysChloris.Radiation: direct_shortwave_trees, direct_shortwave_surfaces
using UrbanTethysChloris.ModelComponents.Parameters:
    HeightDependentVegetationParameters, initialize_heightdependent_vegetationparameters

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

@testset "direct_shortwave_surfaces" begin
    tree_data = Dict{String,Any}(
        "LAI" => 2.0,
        "SAI" => 0.5,
        "hc" => 1.0,
        "h_disp" => 0.7,
        "d_leaf" => 0.02,
        "CASE_ROOT" => 1,
        "ZR95" => 0.5,
        "ZR50" => 0.3,
        "ZRmax" => 1.0,
        "Rrootl" => 5e-4,
        "PsiL50" => -2.0,
        "PsiX50" => -4.0,
        "FI" => 0.08,
        "Do" => 1000.0,
        "a1" => 10.0,
        "go" => 0.01,
        "CT" => 3,
        "DSE" => 0.649,
        "Ha" => 72000.0,
        "gmes" => 0.02,
        "rjv" => 0.01,
        "Kopt" => 0.5,
        "Knit" => 0.5,
        "Vmax" => 50.0,
        "mSl" => 5.0,
        "e_rel" => 0.0001,
        "e_relN" => 0.3,
        "Psi_sto_00" => -0.5,
        "Psi_sto_50" => -2.0,
        "Sl" => 10.0,
        "SPARTREE" => 1,
    )

    tree_vegetation = initialize_heightdependent_vegetationparameters(FT, tree_data)

    @testset "MATLAB - trees" begin
        # Both trees are completely shaded
        d_tree = FT(0.136505190311419)
        h_can = FT(1.124567474048443)
        h_tree = FT(0.815051903114187)
        LAIt = FT(5.0)
        trees = true
        r_tree = FT(0.05)
        SWR_dir = FT(1.027251526757701)
        theta_n = FT(-0.593239662921201)
        theta_Z = FT(1.356513615208954)
        w_can = FT(1.0)

        SWRdir_g, SWRdir_wsun, SWRdir_wshd, SWRdir_t = direct_shortwave_surfaces(
            h_can,
            w_can,
            d_tree,
            h_tree,
            r_tree,
            theta_Z,
            theta_n,
            SWR_dir,
            LAIt,
            trees,
            tree_vegetation,
        )

        @test SWRdir_g == 0.0
        @test SWRdir_wsun ≈ 0.740553148140844
        @test SWRdir_wshd == 0.0
        @test SWRdir_t ≈ 0.309476060386343
    end

    @testset "MATLAB - no trees" begin
        # Both trees are completely shaded
        d_tree = NaN
        h_can = FT(1.124567474048443)
        h_tree = NaN
        LAIt = NaN
        trees = false
        r_tree = NaN
        SWR_dir = FT(1.832071973971378)
        theta_n = FT(-0.109441709557626)
        theta_Z = FT(1.234325341923081)
        w_can = FT(1.0)

        SWRdir_g, SWRdir_wsun, SWRdir_wshd, SWRdir_t = direct_shortwave_surfaces(
            h_can,
            w_can,
            d_tree,
            h_tree,
            r_tree,
            theta_Z,
            theta_n,
            SWR_dir,
            LAIt,
            trees,
            tree_vegetation,
        )

        @test SWRdir_g ≈ 1.188703834894062
        @test SWRdir_wsun ≈ 0.572102745210290
        @test SWRdir_wshd == 0.0
        @test SWRdir_t == 0.0
    end

    @testset "MATLAB - energy conservation check" begin
        # Configuration that causes initial energy imbalance
        d_tree = FT(0.2)
        h_can = FT(2.0)  # Taller canyon
        h_tree = FT(1.5) # Taller trees
        LAIt = FT(0.5)   # Lower LAI for more direct radiation penetration
        trees = true
        r_tree = FT(0.3) # Larger tree radius
        SWR_dir = FT(1000.0)
        theta_n = FT(π/4) # 45 degree angle
        theta_Z = FT(π/6) # 30 degree zenith
        w_can = FT(1.0)

        SWRdir_g, SWRdir_wsun, SWRdir_wshd, SWRdir_t = direct_shortwave_surfaces(
            h_can,
            w_can,
            d_tree,
            h_tree,
            r_tree,
            theta_Z,
            theta_n,
            SWR_dir,
            LAIt,
            trees,
            tree_vegetation,
        )

        @test SWRdir_g ≈ 142.9126064697672
        @test SWRdir_wsun ≈ 326.801691344312
        @test SWRdir_wshd == 0.0
        @test SWRdir_t ≈ 53.975810276011444
    end
end
