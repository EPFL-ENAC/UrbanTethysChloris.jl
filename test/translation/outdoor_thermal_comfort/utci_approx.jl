using Test
using MAT
using UrbanTethysChloris.OutdoorThermalComfort: utci_approx
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("OTC.UTCI_approx.json")

@testset "MATLAB" begin
    UTCI = utci_approx(
        input_vars["Ta"], input_vars["RH"], input_vars["Tmrt"], input_vars["va"]
    )

    @test UTCI ≈ output_vars["UTCI_approx"]
end
