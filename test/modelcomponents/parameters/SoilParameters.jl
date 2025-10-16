using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    VegetatedSoilParameters,
    SoilParameters,
    initialize_vegetated_soilparameters,
    initialize_wall_soilparameters,
    initialize_soil_parameters,
    VegetationParameters

using ....TestUtils: create_height_dependent_vegetation_parameters

FT = Float64

vegground = create_height_dependent_vegetation_parameters(
    FT; CASE_ROOT=1, ZR95=[250.0], ZR50=[NaN], ZRmax=[NaN]
)

vegroof = create_height_dependent_vegetation_parameters(
    FT; CASE_ROOT=1, ZR95=[70.0], ZR50=[NaN], ZRmax=[NaN]
)

vegtree = create_height_dependent_vegetation_parameters(
    FT; CASE_ROOT=1, ZR95=[1000.0], ZR50=[NaN], ZRmax=[NaN]
)

vegetation_params = VegetationParameters{FT}(; roof=vegroof, ground=vegground, tree=vegtree)

soil_params = Dict{String,Any}()

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
    "dz1" => 0.1,
    "dz2" => 0.1,
    "Zs" => [0.0, 10.0, 20.0, 50.0, 100.0],
    "FixSM" => true,
    "FixSM_LayerStart" => 1,
    "FixSM_LayerEnd" => 4,
)

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
    "Zs" => [
        0.0,
        10.0,
        20.0,
        50.0,
        100.0,
        150.0,
        200.0,
        300.0,
        400.0,
        600.0,
        800.0,
        1000.0,
        1500.0,
        2000.0,
    ],
    "FixSM" => true,
    "FixSM_LayerStart" => 6,
    "FixSM_LayerEnd" => 13,
)

# Wall soil parameters
soil_params["wall"] = Dict{String,Any}("dz1" => 0.1, "dz2" => 0.1)

# Tree interception parameter
soil_params["Sp_In_T"] = 0.2

@testset "VegetatedSoilProperties initialization" begin
    # Test roof soil parameters
    roof_params = soil_params["roof"]
    roof_parameters = initialize_vegetated_soilparameters(FT, roof_params, vegroof, vegroof)

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
    @test roof_parameters.dz1 == FT(0.1)
    @test roof_parameters.dz2 == FT(0.1)
    @test roof_parameters.Zs == [FT(0.0), FT(10.0), FT(20.0), FT(50.0), FT(100.0)]
    @test roof_parameters.ms == 4
    @test roof_parameters.FixSM == true
    @test roof_parameters.FixSM_LayerStart == 1
    @test roof_parameters.FixSM_LayerEnd == 4

    # Test ground soil parameters
    ground_params = soil_params["ground"]
    ground_parameters = initialize_vegetated_soilparameters(
        FT, ground_params, vegtree, vegground
    )

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
    @test ground_parameters.Zs == [
        FT(0.0),
        FT(10.0),
        FT(20.0),
        FT(50.0),
        FT(100.0),
        FT(150.0),
        FT(200.0),
        FT(300.0),
        FT(400.0),
        FT(600.0),
        FT(800.0),
        FT(1000.0),
        FT(1500.0),
        FT(2000.0),
    ]
    @test ground_parameters.ms == 13
    @test ground_parameters.FixSM == true
    @test ground_parameters.FixSM_LayerStart == 6
end

@testset "Extraneous VegetatedSoilParameters field" begin
    extra_params = copy(soil_params["ground"])
    extra_params["extra_field"] = 0.5
    @test_throws ArgumentError initialize_vegetated_soilparameters(
        FT, extra_params, vegroof, vegroof
    )
end

@testset "WallSoilProperties initialization" begin
    # Test wall soil parameters
    wall_params = soil_params["wall"]
    wall_parameters = initialize_wall_soilparameters(FT, wall_params)

    @test wall_parameters.dz1 == FT(0.1)
    @test wall_parameters.dz2 == FT(0.1)
end

@testset "SoilProperties initialization" begin
    @test_nowarn initialize_soil_parameters(FT, soil_params, vegetation_params)

    soil_parameters = initialize_soil_parameters(FT, soil_params, vegetation_params)
    @test soil_parameters.Sp_In_T == FT(0.2)
end
