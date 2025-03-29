using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    LocationSpecificThermalProperties,
    ThermalProperties,
    initialize_thermalproperties,
    initialize_locationspecific_thermalproperties

FT = Float64

@testset "LocationSpecificThermalProperties" begin
    params = Dict{String,Any}("lan_dry" => 0.67, "cv_s" => 1.0e6)
    prop = initialize_locationspecific_thermalproperties(FT, params)

    @test prop.lan_dry == FT(0.67)
    @test prop.cv_s == FT(1.0e6)
end

@testset "ThermalProperties initialization" begin
    params = Dict{String,Any}(
        "roof" => Dict{String,Any}("lan_dry" => 0.67, "cv_s" => 1.0e6),
        "ground" => Dict{String,Any}("lan_dry" => 1.2, "cv_s" => 1.5e6),
        "wall" => Dict{String,Any}("lan_dry" => 0.67, "cv_s" => 1.0e6),
    )
    prop = initialize_thermalproperties(FT, params)

    @test prop.roof.lan_dry == FT(0.67)
    @test prop.ground.lan_dry == FT(1.2)
    @test prop.wall.lan_dry == FT(0.67)
end

@testset "Invalid keys" begin
    params = Dict{String,Any}(
        "roof" => Dict{String,Any}("lan_dry" => 0.67, "cv_s" => 1.0e6),
        "invalid" => Dict{String,Any}("lan_dry" => 1.2, "cv_s" => 1.5e6),
    )
    @test_throws ArgumentError initialize_thermalproperties(FT, params)
end
