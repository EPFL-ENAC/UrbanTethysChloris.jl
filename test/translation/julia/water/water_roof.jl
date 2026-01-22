using Test
using MAT
using UrbanTethysChloris.Water: water_roof
using UrbanTethysChloris.ModelComponents.Parameters
using UrbanTethysChloris.ModelComponents.ForcingInputs
using UrbanTethysChloris.ModelComponents.ModelVariables
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("water_functions.WaterRoof.json")

Eroof_imp = input_vars["Eroof_imp"]
Eroof_veg = FT(input_vars["Eroof_veg"])
Eroof_ground = FT(input_vars["Eroof_ground"])
Eroof_soil_in = input_vars["Eroof_soil"]
TEroof_veg_in = FT(input_vars["TEroof_veg"])

FractionsRoof = LocationSpecificSurfaceFractions(FT, input_vars["FractionsRoof"])
ParSoilRoof = VegetatedSoilParameters(FT, input_vars["ParSoilRoof"], "_R")
ParVegRoof = HeightDependentVegetationParameters(FT, input_vars["ParVegRoof"])

ParCalculation = (
    dth=FT(input_vars["ParCalculation"]["dth"]), row=FT(input_vars["ParCalculation"]["row"])
)

@testset "Model variables" begin
    MeteoData = MeteorologicalInputs(FT, input_vars["MeteoData"])
    Int_ittm = Interception(FT, input_vars["Int_ittm"])
    Owater_ittm = Owater(FT, input_vars["Owater_ittm"])
    Runon_ittm = Runon(FT, input_vars["Runon_ittm"])
    Anthropogenic = AnthropogenicInputs(FT, input_vars["Anthropogenic"])

    q_runon_imp, In_imp, dIn_imp_dt, Lk_imp, q_runon_veg, In_veg, dIn_veg_dt, q_runon_ground, In_ground, dIn_ground_dt, dV_dt, f_ground, V, O, OS, Lk, Psi_s, Exwat, Rd, TEroof_veg_out, Eroof_soil_out, runoff, runon, WBalance_In_imp, WBalance_In_veg, WBalance_In_ground, WBalance_soil, WBalance_imp_tot, WBalance_veg_tot, WBalance_tot = water_roof(
        Eroof_imp,
        Eroof_veg,
        Eroof_ground,
        Eroof_soil_in,
        TEroof_veg_in,
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        FractionsRoof,
        ParSoilRoof,
        ParCalculation,
        ParVegRoof,
        Anthropogenic,
    )

    @test q_runon_imp ≈ output_vars["q_runon_imp"]
    @test In_imp ≈ output_vars["In_imp"]
    @test dIn_imp_dt ≈ output_vars["dIn_imp_dt"]
    @test Lk_imp ≈ output_vars["Lk_imp"]
    @test q_runon_veg ≈ output_vars["q_runon_veg"]
    @test In_veg ≈ output_vars["In_veg"]
    @test dIn_veg_dt ≈ output_vars["dIn_veg_dt"]
    @test q_runon_ground ≈ output_vars["q_runon_ground"]
    @test In_ground ≈ output_vars["In_ground"]
    @test dIn_ground_dt ≈ output_vars["dIn_ground_dt"]
    @test dV_dt ≈ output_vars["dV_dt"]
    @test f_ground ≈ output_vars["f_ground"]
    @test V ≈ vec(output_vars["V"]) atol = 0.05
    @test O ≈ vec(output_vars["O"]) atol = 0.05
    @test OS ≈ output_vars["OS"] atol = 0.05
    @test Lk ≈ output_vars["Lk"] atol = 0.05
    @test Psi_s ≈ output_vars["Psi_s"] atol = 0.05
    # @test all(Exwat .≈ vec(output_vars["Exwat"]))
    @test Rd == output_vars["Rd"]
    @test TEroof_veg_out == output_vars["TEroof_veg"]
    @test Eroof_soil_out ≈ output_vars["Eroof_soil"]
    @test runoff ≈ output_vars["Runoff"]
    @test runon ≈ output_vars["Runon_ittm"]
    @test WBalance_In_imp == output_vars["WBalance_In_imp"]
    @test WBalance_In_veg == output_vars["WBalance_In_veg"]
    @test WBalance_In_ground == output_vars["WBalance_In_ground"]
    @test WBalance_soil ≈ output_vars["WBalance_soil"] atol=1e-14
    @test WBalance_imp_tot == output_vars["WBalance_imp_tot"]
    @test WBalance_veg_tot ≈ output_vars["WBalance_veg_tot"] atol=1e-14
    @test WBalance_tot ≈ output_vars["WBalance_tot"] atol=1e-14
end

@testset "Named Tuples" begin
    MeteoData = (Rain=FT(input_vars["MeteoData"]["Rain"]),)
    Int_ittm = (
        IntRoofImp=FT(input_vars["Int_ittm"]["IntRoofImp"]),
        IntRoofVegPlant=FT(input_vars["Int_ittm"]["IntRoofVegPlant"]),
        IntRoofVegGround=FT(input_vars["Int_ittm"]["IntRoofVegGround"]),
    )
    Owater_ittm = (OwRoofSoilVeg=vec(input_vars["Owater_ittm"]["OwRoofSoilVeg"]),)
    Runon_ittm = (RunonRoofTot=FT(input_vars["Runon_ittm"]["RunonRoofTot"]),)

    Anthropogenic = (Waterf_roof=input_vars["Anthropogenic"]["Waterf_roof"],)

    q_runon_imp, In_imp, dIn_imp_dt, Lk_imp, q_runon_veg, In_veg, dIn_veg_dt, q_runon_ground, In_ground, dIn_ground_dt, dV_dt, f_ground, V, O, OS, Lk, Psi_s, Exwat, Rd, TEroof_veg_out, Eroof_soil_out, runoff, runon, WBalance_In_imp, WBalance_In_veg, WBalance_In_ground, WBalance_soil, WBalance_imp_tot, WBalance_veg_tot, WBalance_tot = water_roof(
        Eroof_imp,
        Eroof_veg,
        Eroof_ground,
        Eroof_soil_in,
        TEroof_veg_in,
        MeteoData,
        Int_ittm,
        Owater_ittm,
        Runon_ittm,
        FractionsRoof,
        ParSoilRoof,
        ParCalculation,
        ParVegRoof,
        Anthropogenic,
    )

    @test q_runon_imp ≈ output_vars["q_runon_imp"]
    @test In_imp ≈ output_vars["In_imp"]
    @test dIn_imp_dt ≈ output_vars["dIn_imp_dt"]
    @test Lk_imp ≈ output_vars["Lk_imp"]
    @test q_runon_veg ≈ output_vars["q_runon_veg"]
    @test In_veg ≈ output_vars["In_veg"]
    @test dIn_veg_dt ≈ output_vars["dIn_veg_dt"]
    @test q_runon_ground ≈ output_vars["q_runon_ground"]
    @test In_ground ≈ output_vars["In_ground"]
    @test dIn_ground_dt ≈ output_vars["dIn_ground_dt"]
    @test dV_dt ≈ output_vars["dV_dt"]
    @test f_ground ≈ output_vars["f_ground"]
    @test V ≈ vec(output_vars["V"]) atol = 0.05
    @test O ≈ vec(output_vars["O"]) atol = 0.05
    @test OS ≈ output_vars["OS"] atol = 0.05
    @test Lk ≈ output_vars["Lk"] atol = 0.05
    @test Psi_s ≈ output_vars["Psi_s"] atol = 0.05
    # @test all(Exwat .≈ vec(output_vars["Exwat"]))
    @test Rd == output_vars["Rd"]
    @test TEroof_veg_out == output_vars["TEroof_veg"]
    @test Eroof_soil_out ≈ output_vars["Eroof_soil"]
    @test runoff ≈ output_vars["Runoff"]
    @test runon ≈ output_vars["Runon_ittm"]
    @test WBalance_In_imp == output_vars["WBalance_In_imp"]
    @test WBalance_In_veg == output_vars["WBalance_In_veg"]
    @test WBalance_In_ground == output_vars["WBalance_In_ground"]
    @test WBalance_soil ≈ output_vars["WBalance_soil"] atol=1e-14
    @test WBalance_imp_tot == output_vars["WBalance_imp_tot"]
    @test WBalance_veg_tot ≈ output_vars["WBalance_veg_tot"] atol=1e-14
    @test WBalance_tot ≈ output_vars["WBalance_tot"] atol=1e-14
end
