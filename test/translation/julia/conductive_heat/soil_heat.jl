using Test
using MAT
using UrbanTethysChloris.ConductiveHeat: soil_heat

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "conductive_heat_functions.Soil_Heat.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

@testset "MATLAB" begin
    G, Tdp = soil_heat(
        Int(input_vars["dt"]),
        input_vars["Ts"],
        input_vars["Tstm1"],
        input_vars["Tdptm1"],
        input_vars["CTt"],
    )

    @test G ≈ output_vars["G"]
    @test Tdp ≈ output_vars["Tdp"]
end
