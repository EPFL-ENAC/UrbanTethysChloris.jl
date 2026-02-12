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

@testset "MATLAB" begin
    a1 = 5.0
    Cc = 4.420097675643373e+02
    Csl = 400.0
    CT = 4
    Do = 2000.0
    Ds = 1.139931698763772e+02
    DS = 0.649
    FI = 0.040
    gmes = Inf
    go = 0.01
    Ha = 72.0
    IPAR = 8.316346043338774
    Oa = 210000.0
    Pre = 9.954590600000000e+02
    Psi_L = -0.034533250966303
    Psi_sto_00 = -0.5
    Psi_sto_50 = -1.6
    Ts = 12.564639365338962
    Vmax = 40.900673280414473
    ra = 1.992126748077743e+02
    rb = 52.861677231827770
    rjv = 2.1

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

    @test CcF ≈ 346.453047878855
    @test An ≈ 1.001323767012722
    @test rs ≈ 1201.607221875981
    @test Rdark ≈ 0.431845909261024
    @test F755nm ≈ 0.336136823355323
    @test GAM ≈ 2.190225452077582
    @test gsCO2 ≈ 21277.52485967636
end
