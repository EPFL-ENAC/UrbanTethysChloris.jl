using Test
using UrbanTethysChloris.Resistance: photosynthesis_biochemical

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

    # Call the function with the test values
    CcF, An, rs, Rdark, F755nm, GAM, gsCO2 = photosynthesis_biochemical(
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

    @test CcF ≈ 425.579109419476
    @test An ≈ -0.245743457481728
    @test rs ≈ 2512.893583003947
    @test Rdark ≈ 0.245743457481728
    @test isnan(F755nm)
    @test GAM ≈ 1.268521967198623
    @test gsCO2 ≈ 10000
end

@testset "TethysChloris.jl issue #136" begin
    CT = 3
    Cc = 200.0
    Csl = 338.43
    DS = 0.656
    Do = 1000.0
    Ds = 310.5037657846657
    FI = 0.081
    Ha = 55.0
    IPAR = 1.602272776775221
    Oa = 210000.0
    Pre = 948.7232105699336
    Psi_L = -0.019349689699981
    Psi_sto_00 = -0.5
    Psi_sto_50 = -3.0
    Ta = 4.1
    Ts = 4.1
    Vmax = 83.571947364501554
    a1 = 7.0
    gmes = Inf
    go = 0.01
    ra = 20.725120860990302
    rb = 15.231030461445188
    rjv = 2.4

    CcF, An, rs, Rdark, F755nm, GAM, gsCO2 = photosynthesis_biochemical(
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

    @test CcF ≈ 329.4521734704098
end
