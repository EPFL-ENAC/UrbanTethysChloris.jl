using Test
using MAT
using UrbanTethysChloris.ConductiveHeat: soil_heat
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("conductive_heat_functions.Soil_Heat.json")

@testset "MATLAB" begin
    G, Tdp = soil_heat(
        Int(input_vars["dt"]),
        input_vars["Ts"],
        input_vars["Tstm1"],
        FT(input_vars["Tdptm1"]),
        input_vars["CTt"],
    )

    @test G ≈ output_vars["G"]
    @test Tdp ≈ output_vars["Tdp"]
end
