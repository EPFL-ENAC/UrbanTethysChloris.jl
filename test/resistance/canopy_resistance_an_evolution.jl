using Test
using UrbanTethysChloris.Resistance: canopy_resistance_an_evolution

FT = Float64

@testset "MATLAB" begin
    a1 = 6.0
    Ca = 400.0
    Citm1_shd = 400.0
    Citm1_sun = 400.0
    CT = 3
    Do = 1000.0
    Ds = 43.591516900888905
    DS = 0.656
    e_rel = 1.0
    e_relN = 1.0
    FI = 0.081
    Fshd = 0.429203837488152
    Fsun = 0.570796162511848
    gmes = Inf
    go = 0.01
    Ha = 55.0
    Knit = 0.15
    Kopt = 0.5
    LAI = 2.5
    mSl = 0.0
    Oa = 210000.0
    Opt_CR = eps(FT)
    PAR_shd = 0.0
    PAR_sun = 0.0
    Pre = 948.7232105699336
    Psi_L = 0.0
    Psi_sto_50 = -3.0
    Psi_sto_99 = -0.5
    ra = 100.0
    rb = 50.0
    rjv = 2.4
    Sl = 0.035
    Ta = 3.899999999999977
    Ts = 3.899999999999977
    Vmax = 68.0

    # Call the function
    rs_sun, rs_shd, Ci_sun, Ci_shd, An, Rdark, Lpho, SIF, DCi = canopy_resistance_an_evolution(
        PAR_sun,
        PAR_shd,
        LAI,
        Kopt,
        Knit,
        Fsun,
        Fshd,
        Citm1_sun,
        Citm1_shd,
        Ca,
        ra,
        rb,
        Ts,
        Ta,
        Pre,
        Ds,
        Psi_L,
        Psi_sto_50,
        Psi_sto_99,
        CT,
        Vmax,
        DS,
        Ha,
        FI,
        Oa,
        Do,
        a1,
        go,
        e_rel,
        e_relN,
        gmes,
        rjv,
        mSl,
        Sl,
        Opt_CR,
    )

    @test rs_sun ≈ 2512.893583003947
    @test rs_shd ≈ 2512.893583003947
    @test Ci_sun ≈ 425.579109419476
    @test Ci_shd ≈ 423.3814758259464
    @test An ≈ -0.591704093976337
    @test Rdark ≈ 0.591704093976336
    @test Lpho ≈ -3.644862189844389e-16 atol = 1e-15
    @test isnan(SIF)
    @test DCi ≈ 0 atol = 1e-13
end

@testset "MATLAB" begin
    a1 = 5.0
    Ca = 400.0
    Citm1_shd = 4.343151191623015e+02
    Citm1_sun = 4.420097675643373e+02
    CT = 4
    Do = 2000.0
    Ds = 1.139931698763772e+02
    DS = 0.649
    e_rel = 1.0
    e_relN = 1.0
    FI = 0.040
    Fshd = 0.429203837488152
    Fsun = 0.570796162511848
    gmes = Inf
    go = 0.01
    Ha = 72.0
    Knit = 0.3
    Kopt = 0.5
    LAI = 2.5
    mSl = 0.0
    Oa = 210000.0
    Opt_CR = eps(FT)
    PAR_shd = 8.923519089201031
    PAR_sun = 11.867346019145904
    Pre = 9.954590600000000e+02
    Psi_L = -0.034533250966303
    Psi_sto_50 = -1.6
    Psi_sto_99 = -0.5
    ra = 1.992126748077743e+02
    rb = 52.861677231827770
    rjv = 2.1
    Sl = 0.025
    Ta = 7.495203358980120
    Ts = 12.564639365338962
    Vmax = 54.0

    # Call the function
    rs_sun, rs_shd, Ci_sun, Ci_shd, An, Rdark, Lpho, SIF, DCi = canopy_resistance_an_evolution(
        PAR_sun,
        PAR_shd,
        LAI,
        Kopt,
        Knit,
        Fsun,
        Fshd,
        Citm1_sun,
        Citm1_shd,
        Ca,
        ra,
        rb,
        Ts,
        Ta,
        Pre,
        Ds,
        Psi_L,
        Psi_sto_50,
        Psi_sto_99,
        CT,
        Vmax,
        DS,
        Ha,
        FI,
        Oa,
        Do,
        a1,
        go,
        e_rel,
        e_relN,
        gmes,
        rjv,
        mSl,
        Sl,
        Opt_CR,
    )

    @test rs_sun ≈ 1050.643824956689
    @test rs_shd ≈ 1016.860647229441
    @test Ci_sun ≈ 352.3737674918041
    @test Ci_shd ≈ 351.2513345148912
    @test An ≈ 2.559602404607898
    @test Rdark ≈ 1.002774256566948
    @test Lpho ≈ 1.670754654091003
    @test SIF ≈ 0.849402936405788
    @test DCi ≈ 0 atol = 1e-13
end
