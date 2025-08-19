using Test
using MAT
using UrbanTethysChloris.Soil: soil_parameters_total
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("soil_functions.SoilParametersTotal.mat")

@testset "MATLAB" begin
    Zs, dz, ms, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, EvL_Zs, Inf_Zs, RfH_Zs, RfL_Zs, Zinf, Kbot, Slo_pot, Dz, aR, aTop, rsd, lan_dry, lan_s, cv_s = soil_parameters_total(
        input_vars["Pcla"],
        input_vars["Psan"],
        input_vars["Porg"],
        input_vars["Kfc"],
        input_vars["Phy"],
        Int(input_vars["SPAR"]),
        input_vars["Kbot"],
        Int(input_vars["CASE_ROOT_H"]),
        Int(input_vars["CASE_ROOT_L"]),
        [input_vars["ZR95_H"]],
        [input_vars["ZR95_L"]],
        [input_vars["ZR50_H"]],
        [input_vars["ZR50_L"]],
        [input_vars["ZRmax_H"]],
        [input_vars["ZRmax_L"]],
        vec(input_vars["Zs"]),
    )

    # Compare the result with the expected output
    @test all(Zs .≈ vec(output_vars["Zs"]))
    @test all(dz .≈ vec(output_vars["dz"]))
    @test ms == output_vars["ms"]
    @test all(Osat .≈ output_vars["Osat"])
    @test all(Ohy .≈ output_vars["Ohy"])
    @test all(nVG .≈ output_vars["nVG"])
    @test all(alpVG .≈ output_vars["alpVG"])
    @test all(Ks_Zs .≈ output_vars["Ks_Zs"])
    @test all(L .≈ output_vars["L"])
    @test all(Pe .≈ output_vars["Pe"])
    @test all(O33 .≈ output_vars["O33"])
    @test SPAR == output_vars["SPAR"]
    @test all(EvL_Zs .≈ vec(output_vars["EvL_Zs"]))
    @test all(Inf_Zs .≈ vec(output_vars["Inf_Zs"]))
    @test all(RfH_Zs .≈ output_vars["RfH_Zs"])
    @test all(RfL_Zs .≈ output_vars["RfL_Zs"])
    @test Zinf ≈ output_vars["Zinf"]
    @test isnan(Kbot)
    @test all(Slo_pot .≈ output_vars["Slo_pot"])
    @test all(Dz .≈ vec(output_vars["Dz"]))
    @test all(aR .≈ output_vars["aR"])
    @test all(aTop .≈ output_vars["aTop"])
    @test all(rsd .≈ output_vars["rsd"])
    @test all(lan_dry .≈ vec(output_vars["lan_dry"]))
    @test all(lan_s .≈ vec(output_vars["lan_s"]))
    @test all(cv_s .≈ vec(output_vars["cv_s"]))
end
