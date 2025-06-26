using Test
using MAT
using UrbanTethysChloris.OutdoorThermalComfort: utci_approx

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "OTC.UTCI_approx.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

@testset "MATLAB" begin
    UTCI = utci_approx(
        input_vars["Ta"], input_vars["RH"], input_vars["Tmrt"], input_vars["va"]
    )

    @test UTCI ≈ output_vars["UTCI_approx"]
end
