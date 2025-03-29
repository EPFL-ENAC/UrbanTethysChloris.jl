using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    VegetatedSoilParameters,
    SoilParameters,
    initialize_vegetated_soilparameters,
    initialize_soil_parameters

FT = Float64

soil_params = Dict{String,Any}()

# Roof soil parameters
soil_params["roof"] = Dict{String,Any}(
    "Pcla" => 0.20,
    "Psan" => 0.40,
    "Porg" => 0.025,
    "In_max_imp" => 0.25,
    "In_max_ground" => 10.0,
    "Sp_In" => 0.2,
    "Kimp" => 0.0,
    "Kfc" => 0.2,
    "Phy" => 10000.0,
    "SPAR" => 2,
    "Kbot" => NaN,
)

# Ground soil parameters
soil_params["ground"] = Dict{String,Any}(
    "Pcla" => 0.20,
    "Psan" => 0.40,
    "Porg" => 0.025,
    "In_max_imp" => 0.5,
    "In_max_underveg" => 10.0,
    "In_max_bare" => 10.0,
    "Sp_In" => 0.2,
    "Kimp" => 0.001,
    "Kfc" => 0.2,
    "Phy" => 10000.0,
    "SPAR" => 2,
    "Kbot" => NaN,
)

# Tree interception parameter
soil_params["Sp_In_T"] = 0.2

@testset "VegetatedSoilProperties initialization" begin
    # Test roof soil parameters
    roof_params = soil_params["roof"]
    roof_parameters = initialize_vegetated_soilparameters(FT, roof_params)

    @test roof_parameters.Pcla == FT(0.20)
    @test roof_parameters.Psan == FT(0.40)
    @test roof_parameters.Porg == FT(0.025)
    @test roof_parameters.In_max_imp == FT(0.25)
    @test roof_parameters.In_max_ground == FT(10.0)
    @test roof_parameters.Sp_In == FT(0.2)
    @test roof_parameters.Kimp == FT(0.0)
    @test roof_parameters.Kfc == FT(0.2)
    @test roof_parameters.Phy == FT(10000.0)
    @test roof_parameters.SPAR == 2
    @test isnan(roof_parameters.Kbot)

    # Test ground soil parameters
    ground_params = soil_params["ground"]
    ground_parameters = initialize_vegetated_soilparameters(FT, ground_params)

    @test ground_parameters.Pcla == FT(0.20)
    @test ground_parameters.Psan == FT(0.40)
    @test ground_parameters.Porg == FT(0.025)
    @test ground_parameters.In_max_imp == FT(0.5)
    @test ground_parameters.In_max_underveg == FT(10.0)
    @test ground_parameters.In_max_bare == FT(10.0)
    @test ground_parameters.Sp_In == FT(0.2)
    @test ground_parameters.Kimp == FT(0.001)
    @test ground_parameters.Kfc == FT(0.2)
    @test ground_parameters.Phy == FT(10000.0)
    @test ground_parameters.SPAR == 2
    @test isnan(ground_parameters.Kbot)
end

@testset "SoilProperties initialization" begin
    @test_nowarn initialize_soil_parameters(FT, soil_params)

    soil_parameters = initialize_soil_parameters(FT, soil_params)
    @test soil_parameters.Sp_In_T == FT(0.2)
end

@testset "SoilProperties validation" begin
    # Test that the validation function throws an error if an extraneous key is present
    extraneous_params = copy(soil_params)
    extraneous_params["extraneous"] = 0.0
    @test_throws ArgumentError initialize_soil_parameters(FT, extraneous_params)
end
