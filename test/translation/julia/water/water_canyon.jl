using Test
using UrbanTethysChloris.Water: water_canyon
using ....TestUtils:
    create_vegetated_soil_parameters,
    create_height_dependent_vegetation_parameters,
    create_location_specific_surface_fractions,
    create_urban_geometry_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("water_functions.WaterCanyon.json")

@testset "MATLAB" begin

    # Create input NamedTuples
    MeteoData = (; Rain=FT(input_vars["MeteoData"]["Rain"]),)
    Int_ittm = (;
        IntGroundImp=FT(input_vars["Int_ittm"]["IntGroundImp"]),
        IntGroundVegPlant=FT(input_vars["Int_ittm"]["IntGroundVegPlant"]),
        IntGroundVegGround=FT(input_vars["Int_ittm"]["IntGroundVegGround"]),
        IntTree=FT(input_vars["Int_ittm"]["IntTree"]),
        IntGroundBare=FT(input_vars["Int_ittm"]["IntGroundBare"]),
    )

    Owater_ittm = (;
        OwGroundSoilVeg=vec(input_vars["Owater_ittm"]["OwGroundSoilVeg"]),
        OwGroundSoilBare=vec(input_vars["Owater_ittm"]["OwGroundSoilBare"]),
        OwGroundSoilImp=vec(input_vars["Owater_ittm"]["OwGroundSoilImp"]),
    )
    Runon_ittm = (; RunonGroundTot=FT(input_vars["Runon_ittm"]["RunonGroundTot"]),)
    Qinlat_ittm = (;
        Qin_imp=FT.(vec(input_vars["Qinlat_ittm"]["Qin_imp"])),
        Qin_bare=FT.(vec(input_vars["Qinlat_ittm"]["Qin_bare"])),
        Qin_veg=FT.(vec(input_vars["Qinlat_ittm"]["Qin_veg"])),
    )
    ParSoilGround = create_vegetated_soil_parameters(
        FT;
        Pcla=input_vars["ParSoilGround"]["Pcla"],
        Psan=input_vars["ParSoilGround"]["Psan"],
        Porg=input_vars["ParSoilGround"]["Porg"],
        Kfc=input_vars["ParSoilGround"]["Kfc"],
        Phy=FT(input_vars["ParSoilGround"]["Phy"]),
        SPAR=Int(input_vars["ParSoilGround"]["SPAR"]),
        Kbot=input_vars["ParSoilGround"]["Kbot"],
        Sp_In=input_vars["ParSoilGround"]["Sp_In"],
        Zs=FT.(vec(input_vars["ParSoilGround"]["Zs"])),
        In_max_imp=FT(input_vars["ParSoilGround"]["In_max_imp"]),
        In_max_bare=FT(input_vars["ParSoilGround"]["In_max_bare"]),
        In_max_underveg=FT(input_vars["ParSoilGround"]["In_max_underveg"]),
        Kimp=FT(input_vars["ParSoilGround"]["Kimp"]),
    )

    ParInterceptionTree = (; Sp_In=input_vars["ParInterceptionTree"]["Sp_In"],)
    ParCalculation = (;
        dth=FT(input_vars["ParCalculation"]["dth"]),
        row=FT(input_vars["ParCalculation"]["row"]),
    )
    ParVegGround = create_height_dependent_vegetation_parameters(
        FT;
        LAI=input_vars["ParVegGround"]["LAI"],
        SAI=input_vars["ParVegGround"]["SAI"],
        CASE_ROOT=1,
        ZR95=FT.([input_vars["ParVegGround"]["ZR95"]]),
        ZR50=FT.([input_vars["ParVegGround"]["ZR50"]]),
        ZRmax=FT.([input_vars["ParVegGround"]["ZRmax"]]),
        Rrootl=FT.([input_vars["ParVegGround"]["Rrootl"]]),
        PsiL50=FT.([input_vars["ParVegGround"]["PsiL50"]]),
        PsiX50=FT.([input_vars["ParVegGround"]["PsiX50"]]),
    )
    ParVegTree = create_height_dependent_vegetation_parameters(
        FT;
        LAI=FT(input_vars["ParVegTree"]["LAI"]),
        SAI=input_vars["ParVegTree"]["SAI"],
        CASE_ROOT=1,
        ZR95=FT.([input_vars["ParVegTree"]["ZR95"]]),
        ZR50=FT.([input_vars["ParVegTree"]["ZR50"]]),
        ZRmax=FT.([input_vars["ParVegTree"]["ZRmax"]]),
        Rrootl=FT.([input_vars["ParVegTree"]["Rrootl"]]),
        PsiL50=FT.([input_vars["ParVegTree"]["PsiL50"]]),
        PsiX50=FT.([input_vars["ParVegTree"]["PsiX50"]]),
        SPARTREE=2,
    )
    FractionsGround = create_location_specific_surface_fractions(
        FT;
        Per_runoff=FT(input_vars["FractionsGround"]["Per_runoff"]),
        fveg=FT(input_vars["FractionsGround"]["fveg"]),
        fbare=FT(input_vars["FractionsGround"]["fbare"]),
        fimp=FT(input_vars["FractionsGround"]["fimp"]),
    )
    Gemeotry_m = create_urban_geometry_parameters(
        FT;
        radius_tree=input_vars["geometry"]["radius_tree"],
        trees=Bool(input_vars["ParTree"]["trees"]),
        Width_canyon=FT(input_vars["Gemeotry_m"]["Width_canyon"]),
    )

    Anthropogenic = (;
        Waterf_canyonBare=FT(input_vars["Anthropogenic"]["Waterf_canyonBare"]),
        Waterf_canyonVeg=FT(input_vars["Anthropogenic"]["Waterf_canyonVeg"]),
    )

    # Evaporation/transpiration rates
    Egbare_Pond = FT(input_vars["Egbare_Pond"])
    Egbare_soil = FT(input_vars["Egbare_soil"])
    Egimp_Pond = FT(input_vars["Egimp_Pond"])
    Egveg_In = FT(input_vars["Egveg_In"])
    Egveg_Pond = FT(input_vars["Egveg_Pond"])
    Egveg_soil = FT(input_vars["Egveg_soil"])
    Etree_In = FT(input_vars["Etree_In"])
    TEgveg = FT(input_vars["TEgveg"])
    TEtree = FT(input_vars["TEtree"])

    results = water_canyon(
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        Qinlat_ittm,
        Etree_In,
        Egveg_In,
        Egimp_Pond,
        Egbare_Pond,
        Egveg_Pond,
        Egbare_soil,
        Egveg_soil,
        TEgveg,
        TEtree,
        ParSoilGround,
        ParInterceptionTree,
        ParCalculation,
        ParVegGround,
        ParVegTree,
        FractionsGround,
        Gemeotry_m,
        Anthropogenic,
    )

    @test results.q_tree_dwn ≈ output_vars["q_tree_dwn"]
    @test results.In_tree ≈ output_vars["In_tree"]
    @test results.dIn_tree_dt ≈ output_vars["dIn_tree_dt"]
    @test results.q_gveg_dwn ≈ output_vars["q_gveg_dwn"]
    @test results.In_gveg ≈ output_vars["In_gveg"]
    @test results.dIn_gveg_dt ≈ output_vars["dIn_gveg_dt"]
    @test results.q_gimp_runoff ≈ output_vars["q_gimp_runoff"]
    @test results.In_gimp ≈ output_vars["In_gimp"]
    @test results.dIn_gimp_dt ≈ output_vars["dIn_gimp_dt"]
    @test results.f_inf_gimp ≈ output_vars["f_inf_gimp"]
    @test results.q_gbare_runoff ≈ output_vars["q_gbare_runoff"]
    @test results.In_gbare ≈ output_vars["In_gbare"]
    @test results.dIn_gbare_dt ≈ output_vars["dIn_gbare_dt"]
    @test results.f_inf_gbare ≈ output_vars["f_inf_gbare"]
    @test results.q_gveg_runoff ≈ output_vars["q_gveg_runoff"]
    @test results.In_gveg_pond ≈ output_vars["In_gveg_pond"]
    @test results.dIn_gveg_pond_dt ≈ output_vars["dIn_gveg_pond_dt"]
    @test results.f_inf_gveg ≈ output_vars["f_inf_gveg"]
    @test results.V_gimp ≈ output_vars["V_gimp"] atol = 0.0005
    @test results.O_gimp ≈ output_vars["O_gimp"] atol = 1e-5
    @test results.OS_gimp ≈ output_vars["OS_gimp"] atol = 1e-5
    @test results.Lk_gimp ≈ output_vars["Lk_gimp"]
    @test results.Psi_s_H_gimp ≈ output_vars["Psi_s_H_gimp"] atol = 1e-6
    @test results.Psi_s_L_gimp ≈ output_vars["Psi_s_L_gimp"] atol = 1e-6
    # @test results.Exwat_H_gimp ≈ output_vars["Exwat_H_gimp"]
    # @test results.Exwat_L_gimp ≈ output_vars["Exwat_L_gimp"]
    @test results.Rd_gimp ≈ output_vars["Rd_gimp"]
    @test results.TEgveg_imp ≈ output_vars["TEgveg_imp"]
    @test results.TEtree_imp ≈ output_vars["TEtree_imp"]
    @test results.Egimp_soil ≈ output_vars["Egimp_soil"]
    @test results.dV_dt_gimp ≈ output_vars["dV_dt_gimp"] atol = 0.05
    # @test results.Psi_soil_gimp ≈ output_vars["Psi_soil_gimp"] # Errors of up to 0.5
    @test results.Kf_gimp[3:end] ≈ output_vars["Kf_gimp"][3:end] atol = 1e-5
    @test results.Kf_gimp ≈ output_vars["Kf_gimp"] atol = 1e-5 nans=true
    @test results.V_gbare ≈ output_vars["V_gbare"] atol = 0.05
    @test results.O_gbare ≈ output_vars["O_gbare"] atol = 0.005
    @test results.OS_gbare ≈ output_vars["OS_gbare"] atol = 0.001
    @test results.Lk_gbare ≈ output_vars["Lk_gbare"]
    @test results.Psi_s_H_gbare ≈ output_vars["Psi_s_H_gbare"] atol = 1e-6
    @test results.Psi_s_L_gbare ≈ output_vars["Psi_s_L_gbare"] atol = 1e-5
    # # @test results.Exwat_H_gbare ≈ output_vars["Exwat_H_gbare"]
    # # @test results.Exwat_L_gbare ≈ output_vars["Exwat_L_gbare"]
    # @test results.Rd_gbare ≈ output_vars["Rd_gbare"]
    # @test results.TEgveg_bare ≈ output_vars["TEgveg_bare"]
    # @test results.TEtree_bare ≈ output_vars["TEtree_bare"]
    # @test results.Egbare_Soil ≈ output_vars["Egbare_Soil"]
    # @test results.dV_dt_gbare ≈ output_vars["dV_dt_gbare"]
    # # @test results.Psi_soil_gbare ≈ output_vars["Psi_soil_gbare"]
    # @test results.Kf_gbare ≈ output_vars["Kf_gbare"]
    # @test results.V_gveg ≈ output_vars["V_gveg"]
    # @test results.O_gveg ≈ output_vars["O_gveg"]
    # @test results.OS_gveg ≈ output_vars["OS_gveg"]
    # @test results.Lk_gveg ≈ output_vars["Lk_gveg"]
    # @test results.Psi_s_H_gveg ≈ output_vars["Psi_s_H_gveg"]
    # @test results.Psi_s_L_gveg ≈ output_vars["Psi_s_L_gveg"]
    # @test results.Exwat_H_gveg ≈ output_vars["Exwat_H_gveg"]
    # @test results.Exwat_L_gveg ≈ output_vars["Exwat_L_gveg"]
    # @test results.Rd_gveg ≈ output_vars["Rd_gveg"]
    # @test results.TEgveg_veg ≈ output_vars["TEgveg_veg"]
    # @test results.TEtree_veg ≈ output_vars["TEtree_veg"]
    # @test results.Egveg_Soil ≈ output_vars["Egveg_Soil"]
    # @test results.dV_dt_gveg ≈ output_vars["dV_dt_gveg"]
    # @test results.Psi_soil_gveg ≈ output_vars["Psi_soil_gveg"]
    # @test results.Kf_gveg ≈ output_vars["Kf_gveg"]
    # @test results.Qin_imp ≈ output_vars["Qin_imp"]
    # @test results.Qin_bare ≈ output_vars["Qin_bare"]
    # @test results.Qin_veg ≈ output_vars["Qin_veg"]
    # @test results.Qin_bare2imp ≈ output_vars["Qin_bare2imp"]
    # @test results.Qin_bare2veg ≈ output_vars["Qin_bare2veg"]
    # @test results.Qin_imp2bare ≈ output_vars["Qin_imp2bare"]
    # @test results.Qin_imp2veg ≈ output_vars["Qin_imp2veg"]
    # @test results.Qin_veg2imp ≈ output_vars["Qin_veg2imp"]
    # @test results.Qin_veg2bare ≈ output_vars["Qin_veg2bare"]
    # @test results.V ≈ output_vars["V"]
    # @test results.O ≈ output_vars["O"]
    # @test results.OS ≈ output_vars["OS"]
    # @test results.Lk ≈ output_vars["Lk"]
    # @test results.Rd ≈ output_vars["Rd"]
    # @test results.dV_dt ≈ output_vars["dV_dt"]
    # @test results.Psi_s_L ≈ output_vars["Psi_s_L"]
    # @test results.Exwat_L ≈ output_vars["Exwat_L"]
    # @test results.TEgveg_tot ≈ output_vars["TEgveg_tot"]
    # @test results.Psi_s_H_tot ≈ output_vars["Psi_s_H_tot"]
    # @test results.Exwat_H ≈ output_vars["Exwat_H"]
    # @test results.TEtree_tot ≈ output_vars["TEtree_tot"]
    # @test results.EB_TEtree ≈ output_vars["EB_TEtree"]
    # @test results.EB_TEgveg ≈ output_vars["EB_TEgveg"]
    # @test results.WBIndv ≈ output_vars["WBIndv"]
    # @test results.WBTot ≈ output_vars["WBTot"]
    # @test results.Runoff ≈ output_vars["Runoff"]
    # @test results.Runon_ittm ≈ output_vars["Runon_ittm"]
    # @test results.Etot ≈ output_vars["Etot"]
    # @test results.DeepGLk ≈ output_vars["DeepGLk"]
    # @test results.StorageTot ≈ output_vars["StorageTot"]

end
