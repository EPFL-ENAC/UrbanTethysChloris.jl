using Test
using UrbanTethysChloris.Soil: soil_parameters_total

FT = Float64

@testset "MATLAB" begin
    # Input parameters matching the MATLAB test case
    Zs = FT.([
        0, 10, 20, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1350, 1500
    ])
    Psan = 0.4
    Pcla = 0.2
    Porg = 0.025
    Kfc = 1.0
    Phy = 10000.0
    SPAR = 2
    Kbot = 0.1
    CASE_ROOT = 1
    ZR95_H = [0.0]
    ZR50_H = [NaN]
    ZR95_L = [250.0]
    ZR50_L = [NaN]
    ZRmax_H = [NaN]
    ZRmax_L = [NaN]

    # Call the function
    outputs = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT,
        CASE_ROOT,
        ZR95_H,
        ZR95_L,
        ZR50_H,
        ZR50_L,
        ZRmax_H,
        ZRmax_L,
        Zs,
    )

    # Extract outputs
    Zs_out,
    dz,
    ms,
    Osat,
    Ohy,
    nVG,
    alpVG,
    Ks_Zs,
    L,
    Pe,
    O33,
    SPAR_out,
    EvL_Zs,
    Inf_Zs,
    RfH_Zs,
    RfL_Zs,
    Zinf,
    Kbot_out,
    Slo_pot,
    Dz,
    aR,
    aTop,
    rsd,
    lan_dry,
    lan_s,
    cv_s = outputs

    # Test key outputs against known MATLAB values
    @test all(Osat .≈ 0.459478244940800)
    @test all(L .≈ 0.186873685240403)
    @test all(Pe .≈ 3.912121076244333)
    @test all(Ks_Zs .≈ 15.475656399419266)
    @test all(O33 .≈ 0.279610164940800)
    @test all(rsd .≈ 1432.382650906880)
    @test all(lan_dry .≈ 0.191773412109454)
    @test all(lan_s .≈ 6.537988142956436)
    @test all(cv_s .≈ 2.226788982559455e6)

    # Test root fractions against known values
    expected_RfH_Zs = zeros(16)
    expected_RfL_Zs = [
        [
            0.119006065336605,
            0.105548911440333,
            0.250280177828380,
            0.260595058073905,
            0.143017600179558,
            0.078489723144784,
            0.043076073375656,
        ]
        zeros(9)
    ]

    @test all(vec(RfH_Zs) .≈ expected_RfH_Zs)
    @test all(vec(RfL_Zs) .≈ expected_RfL_Zs)

    # Test output dimensions and types
    @test length(Zs_out) == length(Zs)
    @test length(dz) == ms
    @test size(RfH_Zs) == (1, ms)
    @test size(RfL_Zs) == (1, ms)
    @test all(x -> isa(x, FT), Osat)
    @test isa(aR, FT)
    @test isa(SPAR_out, Int)
end
