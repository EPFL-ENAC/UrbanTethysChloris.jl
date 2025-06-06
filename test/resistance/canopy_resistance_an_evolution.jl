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
    Opt_CR = 1.0
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
