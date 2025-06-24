using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    initialize_building_energy_model_parameters,
    initialize_indooropticalproperties,
    initialize_thermalbuildingparameters,
    initialize_windowparameters,
    initialize_hvacparameters

data = Dict{String,Any}(
    "indoor_optical" => Dict{String,Any}(
        "abc" => 0.3,
        "abw" => 0.3,
        "abg" => 0.3,
        "abm" => 0.3,
        "ec" => 0.95,
        "eg" => 0.95,
        "ew" => 0.95,
        "em" => 0.95,
    ),
    "thermal" => Dict{String,Any}(
        "IntMassOn" => false,
        "FloorHeight" => 3.0,
        "dzFloor" => 0.2,
        "dzWall" => 0.2,
        "lan_ground_floor" => 1.2,
        "cv_ground_floor" => 1.5e6,
        "lan_floor_IntMass" => 0.67,
        "cv_floor_IntMass" => 1.0e6,
        "lan_wall_IntMass" => 0.67,
        "cv_wall_IntMass" => 1.0e6,
    ),
    "windows" => Dict{String,Any}(
        "WindowsOn" => 1,
        "GlazingRatio" => 0.15,
        "Uvalue" => 4.95,
        "lan_windows" => NaN,
        "cv_glass" => 2.1e6,
        "dztot" => 0.02,
        "SHGC" => 0.8,
        "SolarTransmittance" => 0.6,
        "SolarAbsorptivity" => 0.0,
        "SolarAlbedo" => 0.4,
    ),
    "hvac" => Dict{String,Any}(
        "ACon" => true,
        "Heatingon" => true,
        "TsetpointCooling" => 298.15,
        "TsetpointHeating" => 293.15,
        "RHsetpointCooling" => 60.0,
        "RHsetpointHeating" => NaN,
        "ACH" => 0.5,
        "COPAC" => 3.26,
        "COPHeat" => 0.9,
        "f_ACLatentToQ" => 1.0,
    ),
)

FT = Float64

@testset "IndoorOpticalProperties initialization" begin
    io_input = data["indoor_optical"]
    io_parameters = initialize_indooropticalproperties(FT, io_input)

    @test io_parameters.abc == FT(0.3)
    @test io_parameters.abw == FT(0.3)
    @test io_parameters.abg == FT(0.3)
    @test io_parameters.abm == FT(0.3)
    @test io_parameters.ec == FT(0.95)
    @test io_parameters.eg == FT(0.95)
    @test io_parameters.ew == FT(0.95)
    @test io_parameters.em == FT(0.95)
end

@testset "ThermalBuilding initialization" begin
    tb_input = data["thermal"]
    tb_parameters = initialize_thermalbuildingparameters(FT, tb_input)

    @test tb_parameters.IntMassOn == 0
    @test tb_parameters.FloorHeight ≈ 3.0
    @test tb_parameters.dzFloor ≈ 0.2
    @test tb_parameters.dzWall ≈ 0.2
    @test tb_parameters.lan_ground_floor ≈ 1.2
    @test tb_parameters.cv_ground_floor ≈ 1.5e6
    @test tb_parameters.lan_floor_IntMass ≈ 0.67
    @test tb_parameters.cv_floor_IntMass ≈ 1.0e6
    @test tb_parameters.lan_wall_IntMass ≈ 0.67
    @test tb_parameters.cv_wall_IntMass ≈ 1.0e6
end

@testset "WindowParameters initialization" begin
    wp_input = data["windows"]
    wp_parameters = initialize_windowparameters(FT, wp_input)

    @test wp_parameters.WindowsOn == 1
    @test wp_parameters.GlazingRatio ≈ 0.15
    @test wp_parameters.Uvalue ≈ 4.95
    @test isnan(wp_parameters.lan_windows)
    @test wp_parameters.cv_glass ≈ 2.1e6
    @test wp_parameters.dztot ≈ 0.02
    @test wp_parameters.SHGC ≈ 0.8
    @test wp_parameters.SolarTransmittance ≈ 0.6
    @test wp_parameters.SolarAbsorptivity ≈ 0.0
    @test wp_parameters.SolarAlbedo ≈ 0.4
end

@testset "HVACParameters initialization" begin
    hvac_input = data["hvac"]
    hvac_parameters = initialize_hvacparameters(FT, hvac_input)

    @test hvac_parameters.ACon == true
    @test hvac_parameters.Heatingon == true
    @test hvac_parameters.TsetpointCooling ≈ 298.15
    @test hvac_parameters.TsetpointHeating ≈ 293.15
    @test hvac_parameters.RHsetpointCooling ≈ 60.0
    @test isnan(hvac_parameters.RHsetpointHeating)
    @test hvac_parameters.ACH ≈ 0.5
    @test hvac_parameters.COPAC ≈ 3.26
    @test hvac_parameters.COPHeat ≈ 0.9
    @test hvac_parameters.f_ACLatentToQ ≈ 1.0

    #optional parameters
    @test hvac_parameters.AC_onCool == false
    @test hvac_parameters.AC_onDehum == false
    @test hvac_parameters.MasterOn == false
    @test hvac_parameters.q_RHspCooling == 0
end

@testset "BuildingEnergyModelParameters initialization" begin
    params = initialize_building_energy_model_parameters(FT, data)

    @test_throws ArgumentError initialize_building_energy_model_parameters(
        FT, Dict{String,Any}("indoor_optical" => Dict{String,Any}())
    )

    @test_throws ArgumentError initialize_building_energy_model_parameters(
        FT,
        Dict{String,Any}(
            "indoor_optical" => data["indoor_optical"],
            "thermal" => Dict{String,Any}(),
            "windows" => data["windows"],
            "hvac" => data["hvac"],
        ),
    )
end
