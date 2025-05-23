using Test
using UrbanTethysChloris.Resistance: urban_roughness

FT = Float64

@testset "MATLAB" begin
    # Test inputs
    hc_H = 0.15
    hc_L = 0.0
    Csoil = false
    Croad = false
    Croof = true

    zom, zoh, zom_ground, zoh_ground, disp_h, zom_H, zom_L, zoh_H, zoh_L, d_H, d_L, zom_other = urban_roughness(
        hc_H, hc_L, Csoil, Croad, Croof
    )

    @test zom ≈ 0.01845
    @test zoh ≈ 0.001845
    @test zom_ground ≈ 0.01
    @test zoh_ground ≈ 0.001
    @test disp_h ≈ 0.1
    @test zom_H ≈ 0.01845
    @test zom_L ≈ 0.0
    @test zoh_H ≈ 0.001845
    @test zoh_L ≈ 0.0
    @test d_H ≈ 0.1
    @test d_L ≈ 0.0
    @test zom_other ≈ 0.01
end
