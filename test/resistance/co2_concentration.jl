using Test
using UrbanTethysChloris.Resistance: co2_concentration

FT = Float64

@testset "MATLAB" begin
    a1 = 6.0
    Cc = 400.0
    Csl = 400.0
    CT = 3
    Do = 1000.0
    Ds = 43.591516900888905
    DS = 0.656
    FI = 0.081
    gmes = Inf
    go = 0.01
    Ha = 55.0
    IPAR = 0.0
    Oa = 210000.0
    Pre = 948.7232105699336
    Psi_L = 0.0
    Psi_sto_00 = -0.5
    Psi_sto_50 = -3.0
    ra = 100.0
    rb = 50.0
    rjv = 2.4
    Ts = 3.899999999999977
    Vmax = 58.87593470069573

    DCi = co2_concentration(
        Cc,
        IPAR,
        Csl,
        ra,
        rb,
        Ts,
        Pre,
        Ds,
        Psi_L,
        Psi_sto_50,
        Psi_sto_00,
        CT,
        Vmax,
        DS,
        Ha,
        FI,
        Oa,
        Do,
        a1,
        go,
        gmes,
        rjv,
    )
    @test DCi ≈ -25.579109419476
end
