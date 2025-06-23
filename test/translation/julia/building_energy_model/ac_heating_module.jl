using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: ac_heating_module

FT = Float64
dir = joinpath(@__DIR__, "..", "..", "matlab", "data")
filename = "BuildingEnergyModel.AC_HeatingModule.mat"
input_vars = matread(joinpath(dir, "inputs", filename))
output_vars = matread(joinpath(dir, "outputs", filename))

ParHVAC = (MasterOn=input_vars["ParHVAC"]["MasterOn"],)

@testset "MATLAB" begin
    AC_on, AC_onCool, AC_onDehum, Heat_on, H_AC_Heat, LE_AC_Heat = ac_heating_module(
        Bool(input_vars["AC_on"]),
        Bool(input_vars["Heat_on"]),
        Bool(input_vars["AC_onCool"]),
        Bool(input_vars["AC_onDehum"]),
        ParHVAC,
        input_vars["HbuildIn"],
        input_vars["Hvent"],
        input_vars["Hequip"],
        input_vars["Hpeople"],
        input_vars["dSH_air"],
        input_vars["LEvent"],
        input_vars["LEequip"],
        input_vars["LEpeople"],
        input_vars["dSLE_air"],
    )

    @test AC_on == Bool(output_vars["AC_on"])
    @test AC_onCool == Bool(output_vars["AC_onCool"])
    @test AC_onDehum == Bool(output_vars["AC_onDehum"])
    @test Heat_on == Bool(output_vars["Heat_on"])
    @test H_AC_Heat ≈ output_vars["H_AC_Heat"]
    @test LE_AC_Heat ≈ output_vars["LE_AC_Heat"]
end
