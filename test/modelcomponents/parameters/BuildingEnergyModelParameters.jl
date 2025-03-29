using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    initialize_building_energy_model_parameters,
    initialize_indooropticalproperties,
    initialize_thermalbuildingparameters,
    initialize_windowparameters,
    initialize_hvacparameters

# Test data based on Data_UEHM_site_ZH.m
data = Dict{String,Any}(
    "indoor_optical" => Dict{String,Any}(
        "abc" => 0.3,  # Albedo ceiling
        "abw" => 0.3,  # Albedo wall
        "abg" => 0.3,  # Albedo ground
        "abm" => 0.3,  # Albedo internal mass
        "ec" => 0.95,  # Emissivity ceiling
        "eg" => 0.95,  # Emissivity ground
        "ew" => 0.95,  # Emissivity wall
        "em" => 0.95,  # Emissivity internal mass
    ),
    "thermal" => Dict{String,Any}(
        "IntMassOn" => 0,      # Include building internal mass
        "FloorHeight" => 3.0,   # Average floor height (m)
        "dzFloor" => 0.2,      # Average floor thickness (m)
        "dzWall" => 0.2,       # Average wall thickness (m)
        "lan_ground_floor" => 1.2,      # Ground thermal conductivity (W/m K)
        "cv_ground_floor" => 1.5e6,     # Ground volumetric heat capacity (J/m³ K)
        "lan_floor_IntMass" => 0.67,    # Floor thermal conductivity (W/m K)
        "cv_floor_IntMass" => 1.0e6,    # Floor volumetric heat capacity (J/m³ K)
        "lan_wall_IntMass" => 0.67,     # Wall thermal conductivity (W/m K)
        "cv_wall_IntMass" => 1.0e6,     # Wall volumetric heat capacity (J/m³ K)
    ),
    "windows" => Dict{String,Any}(
        "WindowsOn" => 1,          # Include windows in simulation
        "GlazingRatio" => 0.15,    # Window-to-wall ratio
        "Uvalue" => 4.95,          # U-value (W/m² K)
        "lan_windows" => NaN,      # Thermal conductivity (W/m K)
        "cv_glass" => 2.1e6,       # Volumetric heat capacity (J/m³ K)
        "dztot" => 0.02,          # Total glass thickness (m)
        "SHGC" => 0.8,            # Solar heat gain coefficient
        "SolarTransmittance" => 0.6, # Solar transmittance
        "SolarAbsorptivity" => 0.0,  # Solar absorptivity
        "SolarAlbedo" => 0.4,        # Solar albedo
    ),
    "hvac" => Dict{String,Any}(
        "ACon" => 1,                # Enable AC
        "Heatingon" => 1,           # Enable heating
        "TsetpointCooling" => 298.15, # 25°C
        "TsetpointHeating" => 293.15, # 20°C
        "RHsetpointCooling" => 60.0,  # Cooling RH setpoint (%)
        "RHsetpointHeating" => NaN,   # Heating RH setpoint (%)
        "ACH" => 0.5,               # Air changes per hour
        "COPAC" => 3.26,            # AC coefficient of performance
        "COPHeat" => 0.9,           # Heating coefficient of performance
        "f_ACLatentToQ" => 1.0,     # Fraction of latent heat condensed
    ),
)

FT = Float64

@testset "IndoorOpticalProperties initialization" begin
    # Test roof soil parameters
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

    @test hvac_parameters.ACon == 1
    @test hvac_parameters.Heatingon == 1
    @test hvac_parameters.TsetpointCooling ≈ 298.15
    @test hvac_parameters.TsetpointHeating ≈ 293.15
    @test hvac_parameters.RHsetpointCooling ≈ 60.0
    @test isnan(hvac_parameters.RHsetpointHeating)
    @test hvac_parameters.ACH ≈ 0.5
    @test hvac_parameters.COPAC ≈ 3.26
    @test hvac_parameters.COPHeat ≈ 0.9
    @test hvac_parameters.f_ACLatentToQ ≈ 1.0
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
