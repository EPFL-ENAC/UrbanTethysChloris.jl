using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    SurfaceFractions,
    LocationSpecificSurfaceFractions,
    initialize_surfacefractions,
    initialize_locationspecific_surfacefractions

using UrbanTethysChloris.ModelComponents.Parameters:
    SimpleOpticalProperties,
    VegetatedOpticalProperties,
    initialize_simple_opticalproperties,
    initialize_vegetated_opticalproperties,
    initialize_optical_properties

FT = Float64

optical_params = Dict{String,Any}()
optical_params["wall"] = Dict{String,Any}("albedo" => 0.4, "emissivity" => 0.95)
optical_params["tree"] = Dict{String,Any}("albedo" => 0.2, "emissivity" => 0.95)

LAI_R = 2.5
SAI_R = 0.001
optical_params["roof"] = Dict{String,Any}(
    "aveg" => 0.2, "aimp" => 0.15, "eveg" => 1 - exp(-(LAI_R + SAI_R)), "eimp" => 0.95
)

LAI_G = 2.5
SAI_G = 0.001
optical_params["ground"] = Dict{String,Any}(
    "aveg" => 0.2,
    "abare" => 0.15,
    "aimp" => 0.1,
    "eveg" => 1 - exp(-(LAI_G + SAI_G)),
    "ebare" => 0.95,
    "eimp" => 0.95,
)

fractions_params = Dict{String,Any}()
fractions_params["roof"] = Dict{String,Any}(
    "fveg" => 0.5, "fimp" => 0.5, "Per_runoff" => 1.0
)
fractions_params["ground"] = Dict{String,Any}(
    "fveg" => 0.545, "fbare" => 0.0, "fimp" => 0.455, "Per_runoff" => 0.9
)

fractions = initialize_surfacefractions(FT, fractions_params)

@testset "SimpleOpticalProperties initialization" begin
    wall_params = optical_params["wall"]
    prop = initialize_simple_opticalproperties(FT, wall_params)

    @test prop.albedo == FT(0.4)
    @test prop.emissivity == FT(0.95)
end

@testset "VegetatedOpticalProperties initialization" begin
    # Roof
    roof_params = optical_params["roof"]
    roof_properties = initialize_vegetated_opticalproperties(
        FT, roof_params, fractions.roof
    )

    @test roof_properties.albedo ≈ FT(0.175)
    @test roof_properties.emissivity ≈ FT(0.933998522672952)

    # Ground
    ground_params = optical_params["ground"]
    ground_properties = initialize_vegetated_opticalproperties(
        FT, ground_params, fractions.ground
    )

    @test ground_properties.albedo ≈ FT(0.1545)
    @test ground_properties.emissivity ≈ FT(0.932558389713517)
end

@testset "OpticalProperties initialization" begin
    @test_nowarn initialize_optical_properties(FT, optical_params, fractions)
end
