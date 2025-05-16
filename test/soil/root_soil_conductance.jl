using Test
using UrbanTethysChloris.Soil: root_soil_conductance

FT = Float64

@testset "MATLAB" begin
    dz = 10.0
    Ks = 0.2
    rcyl = 0.002
    Rl = 424.0
    rroot = 0.0005

    Ksr = root_soil_conductance(Ks, Rl, rcyl, rroot, dz)

    @test Ksr ≈ 604605.2412276868
end
