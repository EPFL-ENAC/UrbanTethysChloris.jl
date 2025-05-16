using Test
using UrbanTethysChloris.Soil: soil_water_multilayer

FT = Float64

@testset "MATLAB" begin
    n = 4
    cc = 1
    alpVG = fill(-9.634329919190001e-05, n)
    dz = [10.0, 10.0, 30.0, 50.0]
    EvL_Zs = [1.0, 0.0, 0.0, 0.0]
    Inf_Zs = [1.0, 0.0, 0.0, 0.0]
    Ks_Zs = fill(15.475656399419266, n)
    L = fill(0.186873685240403, n)
    nVG = fill(1.186873685240403, n)
    O33 = fill(0.279610164940800, n)
    Ohy = fill(0.096123264288629, n)
    Osat = fill(0.459478244940800, n)
    Pe = fill(3.912121076244333, n)
    PsiL50_H = fill(0.0, cc)
    PsiL50_L = fill(-4.0, cc)
    PsiX50_H = fill(0.0, cc)
    PsiX50_L = fill(-4.5000, cc)
    RfH_Zs = zeros(FT, cc, n)
    RfL_Zs = zeros(FT, cc, n)
    RfL_Zs[1, :] = [
        0.284979001762539, 0.207810379282956, 0.342622081539846, 0.164602146793880
    ]
    Rrootl_H = fill(0.0, cc)
    Rrootl_L = fill(3000.0, cc)
    SPAR = 2
    V = [-0.961232650478129, -0.961232642886289, -2.883697928658866, -4.806163206839603]
    Zs = [0.0, 10.0, 20.0, 50.0, 100.0]

    O, ZWT, OF, OS, Psi_s_H, Psi_s_L, gsr_H, gsr_L, Exwat_H, Exwat_L, Rd, WTR, POT, OH, OL = soil_water_multilayer(
        V,
        Zs,
        dz,
        n,
        Osat,
        Ohy,
        nVG,
        alpVG,
        Ks_Zs,
        L,
        Pe,
        O33,
        SPAR,
        EvL_Zs,
        Inf_Zs,
        RfH_Zs,
        RfL_Zs,
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    @test all(O .≈ 0.096123264288629)
    @test ZWT ≈ 100
    @test OF ≈ 0.096123264288629
    @test OS ≈ 0.096123264288629
    @test all(isnan.(Psi_s_H))
    @test all(Psi_s_L .≈ -9.999222423640312)
    @test all(gsr_H .≈ 0.0)
    @test all(gsr_L .≈ 0.162155326866212)
    @test all(Exwat_H .≈ 0.0)
    @test all(Exwat_L .≈ 0.0)
    @test Rd ≈ 0.0
    @test all(WTR .≈ 0.0)
    @test all(POT .≈ -1.019367991845057e+06)
    @test all(OH .≈ 0.0)
    @test all(OL .≈ 0.096124572466585)
end
