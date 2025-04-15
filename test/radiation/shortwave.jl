using Test
using UrbanTethysChloris.Radiation:
    direct_shortwave_trees,
    direct_shortwave_surfaces,
    shortwave_absorbed_no_trees,
    shortwave_absorbed_with_trees,
    total_shortwave_absorbed
using UrbanTethysChloris.RayTracing: ViewFactor
using UrbanTethysChloris.ModelComponents.Parameters:
    initialize_heightdependent_vegetationparameters,
    initialize_windowparameters,
    initialize_urbangeometry_parameters,
    initialize_locationspecific_surfacefractions,
    initialize_vegetated_opticalproperties,
    initialize_simple_opticalproperties

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

@testset "shortwave_absorbed_no_trees" begin
    @testset "MATLAB" begin
        agbare = FT(0.2)
        agimp = FT(0.1)
        agveg = FT(0.2)
        awraw = FT(0.3)
        BEM_on = true
        fgbare = FT(0.1)
        fgimp = FT(0.67)
        fgveg = FT(0.23)
        h_can = FT(1.120689655172414)
        SWR_diff = FT(2.009)
        SWR_dir = FT(0.0)
        theta_n = FT(1.279240256283245)
        theta_Z = FT(1.520713609889514)

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

        veg_params = initialize_heightdependent_vegetationparameters(FT, tree_data)

        view_factor = ViewFactor{FT}(
            F_gs_nT=0.381290805156702,
            F_gw_nT=0.309354597421649,
            F_ww_nT=0.447921026139827,
            F_wg_nT=0.276039486930087,
            F_ws_nT=0.276039486930087,
            F_sg_nT=0.381290805156702,
            F_sw_nT=0.309354597421649,
        )

        window_data = Dict{String,Any}(
            "WindowsOn" => 1,
            "GlazingRatio" => 0.15,
            "Uvalue" => 4.95,
            "lan_windows" => NaN,
            "cv_glass" => 2.1e6,
            "dztot" => 0.02,
            "SHGC" => 0.8,
            "SolarTransmittance" => 0.6,
            "SolarAbsorptivity" => 0.0,
            "SolarAlbedo" => 0.4,
        )
        window_params = initialize_windowparameters(FT, window_data)

        w_can = FT(1.0)

        SWRin_nT, SWRout_nT, SWRabs_nT, SWRabsDir_nT, SWRabsDiff_nT, SWREB_nT, albedo_canyon = shortwave_absorbed_no_trees(
            h_can,
            w_can,
            fgveg,
            fgbare,
            fgimp,
            awraw,
            agveg,
            agbare,
            agimp,
            SWR_dir,
            SWR_diff,
            theta_Z,
            theta_n,
            view_factor,
            veg_params,
            window_params,
            BEM_on,
        )
        @testset "SWRin_nT" begin
            @test SWRin_nT.GroundImp ≈ 0.899340759110040
            @test SWRin_nT.GroundBare ≈ 0.899340759110040
            @test SWRin_nT.GroundVeg ≈ 0.899340759110040
            @test SWRin_nT.Tree == 0
            @test SWRin_nT.WallSun ≈ 0.684104921098757
            @test SWRin_nT.WallShade ≈ 0.684104921098757
            @test SWRin_nT.TotalGround ≈ 0.899340759110040
            @test SWRin_nT.TotalCanyon ≈ 2.432679375365874
        end
        @testset "SWRout_nT" begin
            @test SWRout_nT.GroundImp ≈ 0.089934075911004
            @test SWRout_nT.GroundBare ≈ 0.179868151822008
            @test SWRout_nT.GroundVeg ≈ 0.179868151822008
            @test SWRout_nT.Tree == 0
            @test SWRout_nT.WallSun ≈ 0.215493050146108
            @test SWRout_nT.WallShade ≈ 0.215493050146108
            @test SWRout_nT.TotalGround ≈ 0.119612320961635
            @test SWRout_nT.TotalCanyon ≈ 0.602613985082223
        end

        @testset "SWRabs_nT" begin
            @test SWRabs_nT.GroundImp ≈ 0.809406683199036
            @test SWRabs_nT.GroundBare ≈ 0.719472607288032
            @test SWRabs_nT.GroundVeg ≈ 0.719472607288032
            @test SWRabs_nT.Tree == 0
            @test SWRabs_nT.WallSun ≈ 0.468611870952648
            @test SWRabs_nT.WallShade ≈ 0.468611870952648
            @test SWRabs_nT.TotalGround ≈ 0.779728438148405
            @test SWRabs_nT.TotalCanyon ≈ 1.830065390283651
        end

        @testset "SWRabsDir_nT" begin
            @test SWRabsDir_nT.GroundImp == 0.0
            @test SWRabsDir_nT.GroundBare == 0.0
            @test SWRabsDir_nT.GroundVeg == 0.0
            @test SWRabsDir_nT.Tree == 0.0
            @test SWRabsDir_nT.WallSun == 0.0
            @test SWRabsDir_nT.WallShade == 0.0
            @test SWRabsDir_nT.TotalGround == 0.0
            @test SWRabsDir_nT.TotalCanyon == 0.0
        end

        @testset "SWRabsDiff_nT" begin
            @test SWRabsDiff_nT.GroundImp ≈ 0.809406683199036
            @test SWRabsDiff_nT.GroundBare ≈ 0.719472607288032
            @test SWRabsDiff_nT.GroundVeg ≈ 0.719472607288032
            @test SWRabsDiff_nT.Tree == 0
            @test SWRabsDiff_nT.WallSun ≈ 0.468611870952648
            @test SWRabsDiff_nT.WallShade ≈ 0.468611870952648
            @test SWRabsDiff_nT.TotalGround ≈ 0.779728438148405
            @test SWRabsDiff_nT.TotalCanyon ≈ 1.830065390283651
        end

        @test albedo_canyon ≈ 0.089066505583051
    end
end

@testset "shortwave_absorbed_with_trees" begin
    @testset "MATLAB" begin
        agbare = FT(0.2)
        agimp = FT(0.1)
        agveg = FT(0.2)
        at = FT(0.2)
        awraw = FT(0.3)
        BEM_on = true
        d_tree = FT(0.222413793103448)
        fgbare = FT(0.1)
        fgimp = FT(0.67)
        fgveg = FT(0.23)
        h_can = FT(1.120689655172414)
        h_tree = FT(0.812068965517241)
        LAIt = FT(3.0)
        r_tree = FT(0.05)
        SWR_diff = FT(62.041942336163004)
        SWR_dir = FT(41.105057663837009)
        theta_n = FT(1.423974338298374)
        theta_Z = FT(1.326778607312699)

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

        veg_params = initialize_heightdependent_vegetationparameters(FT, tree_data)

        view_factor = ViewFactor{FT}(
            F_gs_T=0.290910447761194,
            F_gt_T=0.102572139303483,
            F_gw_T=0.303258706467662,
            F_ww_T=0.346598392652124,
            F_wt_T=0.147054726368159,
            F_wg_T=0.270600076540375,
            F_ws_T=0.235746804439342,
            F_sg_T=0.290910447761194,
            F_sw_T=0.264199004975124,
            F_st_T=0.180691542288557,
            F_tg_T=0.163248629936597,
            F_tw_T=0.262291660245468,
            F_ts_T=0.287579521301221,
            F_tt_T=0.024588528271246,
        )

        window_data = Dict{String,Any}(
            "WindowsOn" => 1,
            "GlazingRatio" => 0.15,
            "Uvalue" => 4.95,
            "lan_windows" => NaN,
            "cv_glass" => 2.1e6,
            "dztot" => 0.02,
            "SHGC" => 0.8,
            "SolarTransmittance" => 0.6,
            "SolarAbsorptivity" => 0.0,
            "SolarAlbedo" => 0.4,
        )
        window_params = initialize_windowparameters(FT, window_data)

        w_can = FT(1.0)

        SWRin_T, SWRout_T, SWRabs_T, SWRabsDir_T, SWRabsDiff_T, SWREB_T, albedo_canyon = shortwave_absorbed_with_trees(
            h_can,
            w_can,
            h_tree,
            r_tree,
            d_tree,
            fgveg,
            fgbare,
            fgimp,
            awraw,
            agveg,
            agbare,
            agimp,
            at,
            LAIt,
            SWR_dir,
            SWR_diff,
            theta_Z,
            theta_n,
            view_factor,
            veg_params,
            window_params,
            BEM_on,
        )

        @testset "SWRin_T" begin
            @test SWRin_T.GroundImp ≈ 25.988581385433669
            @test SWRin_T.GroundBare ≈ 25.988581385433669
            @test SWRin_T.GroundVeg ≈ 25.988581385433669
            @test SWRin_T.Tree ≈ 24.953548918070030
            @test SWRin_T.WallSun ≈ 55.413414679587369
            @test SWRin_T.WallShade ≈ 22.345374364161280
            @test SWRin_T.TotalGround ≈ 25.988581385433669
            @test SWRin_T.TotalCanyon ≈ 128.8108290578973
        end

        @testset "SWRout_T" begin
            @test SWRout_T.GroundImp ≈ 2.598858138543367
            @test SWRout_T.GroundBare ≈ 5.197716277086734
            @test SWRout_T.GroundVeg ≈ 5.197716277086734
            @test SWRout_T.Tree ≈ 4.990709783614006
            @test SWRout_T.WallSun ≈ 17.455225624070021
            @test SWRout_T.WallShade ≈ 7.038792924710804
            @test SWRout_T.TotalGround ≈ 3.456481324262678
            @test SWRout_T.TotalCanyon ≈ 34.042429963962668
        end

        @testset "SWRabs_T" begin
            @test SWRabs_T.GroundImp ≈ 23.389723246890302
            @test SWRabs_T.GroundBare ≈ 20.790865108346935
            @test SWRabs_T.GroundVeg ≈ 20.790865108346935
            @test SWRabs_T.Tree ≈ 19.962839134456026
            @test SWRabs_T.WallSun ≈ 37.958189055517352
            @test SWRabs_T.WallShade ≈ 15.306581439450476
            @test SWRabs_T.TotalGround ≈ 22.532100061170993
            @test SWRabs_T.TotalCanyon ≈ 94.768399093934619
        end

        @testset "SWRabsDir_T" begin
            @test SWRabsDir_T.GroundImp == 0.0
            @test SWRabsDir_T.GroundBare == 0.0
            @test SWRabsDir_T.GroundVeg == 0.0
            @test SWRabsDir_T.Tree == 0.0
            @test SWRabsDir_T.WallSun ≈ 25.124676015142217
            @test SWRabsDir_T.WallShade == 0.0
            @test SWRabsDir_T.TotalGround == 0.0
            @test SWRabsDir_T.TotalCanyon ≈ 28.156964499728350
        end

        @testset "SWRabsDiff_T" begin
            @test SWRabsDiff_T.GroundImp ≈ 23.389723246890302
            @test SWRabsDiff_T.GroundBare ≈ 20.790865108346935
            @test SWRabsDiff_T.GroundVeg ≈ 20.790865108346935
            @test SWRabsDiff_T.Tree ≈ 19.962839134456026
            @test SWRabsDiff_T.WallSun ≈ 12.833513040375134
            @test SWRabsDiff_T.WallShade ≈ 15.306581439450476
            @test SWRabsDiff_T.TotalGround ≈ 22.532100061170993
            @test SWRabsDiff_T.TotalCanyon ≈ 66.611434594206273
        end

        @test albedo_canyon ≈ 0.081229710084301
    end
end

@testset "total_shortwave_absorbed" begin
    @testset "MATLAB" begin
        BEM_on = true
        SWR_dir = FT(0.362608119924872)
        SWR_diff = FT(32.637391880075121)
        theta_n = FT(-0.808129569557205)
        theta_Z = FT(1.470584035118792)

        geometry_data = Dict{String,Any}(
            "Height_canyon" => 6.5,
            "Width_canyon" => 5.78,
            "Width_roof" => 4.73,
            "Radius_tree" => 0.289,
            "Height_tree" => 4.711,
            "Distance_tree" => 0.789,
            "Hcan_max" => NaN,
            "Hcan_std" => NaN,
            "trees" => true,
            "ftree" => 1.0,
        )
        geometry = initialize_urbangeometry_parameters(FT, geometry_data)

        fractions_ground = initialize_locationspecific_surfacefractions(
            FT,
            Dict{String,Any}(
                "fveg" => 0.545, "fbare" => 0.0, "fimp" => 0.455, "Per_runoff" => 0.9
            ),
        )

        optical_ground = initialize_vegetated_opticalproperties(
            FT,
            Dict{String,Any}(
                "aveg" => 0.2,
                "abare" => 0.15,
                "aimp" => 0.1,
                "eveg" => 1 - exp(-(2.5 + 0.001)),
                "ebare" => 0.95,
                "eimp" => 0.95,
            ),
            fractions_ground,
        )

        optical_wall = initialize_simple_opticalproperties(
            FT, Dict{String,Any}("albedo" => 0.4, "emissivity" => 0.95)
        )

        optical_tree = initialize_simple_opticalproperties(
            FT, Dict{String,Any}("albedo" => 0.2, "emissivity" => 0.994483435579239)
        )

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
        tree_parameters = initialize_heightdependent_vegetationparameters(FT, tree_data)

        view_factor = ViewFactor{FT}(;
            F_gs_nT=0.380308601808867,
            F_gw_nT=0.309845699095567,
            F_ww_nT=0.448951341300808,
            F_wg_nT=0.275524329349596,
            F_ws_nT=0.275524329349596,
            F_sg_nT=0.380308601808867,
            F_sw_nT=0.309845699095567,
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

        window_data = Dict{String,Any}(
            "WindowsOn" => 1,
            "GlazingRatio" => 0.15,
            "Uvalue" => 4.95,
            "lan_windows" => NaN,
            "cv_glass" => 2.1e6,
            "dztot" => 0.02,
            "SHGC" => 0.8,
            "SolarTransmittance" => 0.6,
            "SolarAbsorptivity" => 0.0,
            "SolarAlbedo" => 0.4,
        )
        window_params = initialize_windowparameters(FT, window_data)

        SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t, albedo_canyon = total_shortwave_absorbed(
            geometry,
            SWR_dir,
            SWR_diff,
            theta_n,
            theta_Z,
            fractions_ground,
            optical_ground,
            optical_wall,
            optical_tree,
            tree_parameters,
            view_factor,
            window_params,
            BEM_on,
        )

        @testset "SWRin_t" begin
            @test SWRin_t.GroundImp ≈ 12.985532678378403
            @test SWRin_t.GroundBare == 0
            @test SWRin_t.GroundVeg ≈ 12.985532678378403
            @test SWRin_t.Tree ≈ 11.052578721577198
            @test SWRin_t.WallSun ≈ 10.077773945273417
            @test SWRin_t.WallShade ≈ 9.796053253640592
            @test SWRin_t.TotalGround ≈ 12.985532678378405
            @test SWRin_t.TotalCanyon ≈ 42.279532354122331
        end

        @testset "SWRout_t" begin
            @test SWRout_t.GroundImp ≈ 1.298553267837840
            @test SWRout_t.GroundBare == 0
            @test SWRout_t.GroundVeg ≈ 2.597106535675681
            @test SWRout_t.Tree ≈ 2.210515744315440
            @test SWRout_t.WallSun ≈ 4.031109578109367
            @test SWRout_t.WallShade ≈ 3.918421301456237
            @test SWRout_t.TotalGround ≈ 2.006264798809464
            @test SWRout_t.TotalCanyon ≈ 12.334956664509845
        end

        @testset "SWRabs_t" begin
            @test SWRabs_t.GroundImp ≈ 11.686979410540562
            @test SWRabs_t.GroundBare == 0
            @test SWRabs_t.GroundVeg ≈ 10.388426142702723
            @test SWRabs_t.Tree ≈ 8.842062977261758
            @test SWRabs_t.WallSun ≈ 6.046664367164050
            @test SWRabs_t.WallShade ≈ 5.877631952184355
            @test SWRabs_t.TotalGround ≈ 10.979267879568940
            @test SWRabs_t.TotalCanyon ≈ 29.944575689612485
        end

        @testset "SWRabsDir_t" begin
            @test SWRabsDir_t.GroundImp == 0.0
            @test SWRabsDir_t.GroundBare == 0.0
            @test SWRabsDir_t.GroundVeg == 0.0
            @test SWRabsDir_t.Tree == 0.0
            @test SWRabsDir_t.WallSun ≈ 0.193465378446070
            @test SWRabsDir_t.WallShade == 0.0
            @test SWRabsDir_t.TotalGround == 0.0
            @test SWRabsDir_t.TotalCanyon ≈ 0.217564871954923
        end

        @testset "SWRabsDiff_t" begin
            @test SWRabsDiff_t.GroundImp ≈ 11.686979410540562
            @test SWRabsDiff_t.GroundBare == 0
            @test SWRabsDiff_t.GroundVeg ≈ 10.388426142702723
            @test SWRabsDiff_t.Tree ≈ 8.842062977261758
            @test SWRabsDiff_t.WallSun ≈ 5.853198988717979
            @test SWRabsDiff_t.WallShade ≈ 5.877631952184355
            @test SWRabsDiff_t.TotalGround ≈ 10.979267879568940
            @test SWRabsDiff_t.TotalCanyon ≈ 29.727010817657561
        end

        @test albedo_canyon ≈ 0.092588615466288
    end
end
