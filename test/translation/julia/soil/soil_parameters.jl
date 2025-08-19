using Test
using MAT
using UrbanTethysChloris.Soil: soil_parameters
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("soil_functions.Soil_parameters.json")

@testset "MATLAB" begin
    Osat, L, Pe, Ks, O33, rsd, lan_dry, lan_s, cv_s, K_usle = soil_parameters(
        input_vars["Psan"], input_vars["Pcla"], input_vars["Porg"]
    )

    @test Osat ≈ output_vars["Osat"]
    @test L ≈ output_vars["L"]
    @test Pe ≈ output_vars["Pe"]
    @test Ks ≈ output_vars["Ks"]
    @test O33 ≈ output_vars["O33"]
    @test rsd ≈ output_vars["rsd"]
    @test lan_dry ≈ output_vars["lan_dry"]
    @test lan_s ≈ output_vars["lan_s"]
    @test cv_s ≈ output_vars["cv_s"]
    @test K_usle ≈ output_vars["K_usle"]

    # Test error conditions
    @test_throws ArgumentError soil_parameters(-0.1, 0.5, 0.1)
    @test_throws ArgumentError soil_parameters(0.5, 1.2, 0.1)
    @test_throws ArgumentError soil_parameters(0.5, 0.3, 0.3)
end
