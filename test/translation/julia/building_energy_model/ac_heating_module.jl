using Test
using MAT
using UrbanTethysChloris.BuildingEnergyModel: ac_heating_module
using ....TestUtils: create_hvac_parameters, load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("BuildingEnergyModel.AC_HeatingModule.json")

ParHVAC = create_hvac_parameters(FT; MasterOn=Bool(input_vars["ParHVAC"]["MasterOn"]))

@testset "MATLAB" begin
    AC_on, AC_onCool, AC_onDehum, Heat_on, H_AC_Heat, LE_AC_Heat = ac_heating_module(
        Bool(input_vars["AC_on"]),
        Bool(input_vars["Heat_on"]),
        Bool(input_vars["AC_onCool"]),
        Bool(input_vars["AC_onDehum"]),
        ParHVAC,
        FT(input_vars["HbuildIn"]),
        FT(input_vars["Hvent"]),
        FT(input_vars["Hequip"]),
        FT(input_vars["Hpeople"]),
        FT(input_vars["dSH_air"]),
        FT(input_vars["LEvent"]),
        FT(input_vars["LEequip"]),
        FT(input_vars["LEpeople"]),
        FT(input_vars["dSLE_air"]),
    )

    @test AC_on == Bool(output_vars["AC_on"])
    @test AC_onCool == Bool(output_vars["AC_onCool"])
    @test AC_onDehum == Bool(output_vars["AC_onDehum"])
    @test Heat_on == Bool(output_vars["Heat_on"])
    @test H_AC_Heat ≈ output_vars["H_AC_Heat"]
    @test LE_AC_Heat ≈ output_vars["LE_AC_Heat"]
end
