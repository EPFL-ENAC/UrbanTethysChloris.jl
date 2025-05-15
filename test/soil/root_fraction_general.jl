using Test
using UrbanTethysChloris.Soil: root_fraction_general

FT = Float64

@testset "MATLAB" begin
    Zs = FT.([
        0, 10, 20, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1350, 1500
    ])

    CASE_ROOT = 1
    ZR95_H = [0.0]
    ZR50_H = [NaN]
    ZR95_L = [250.0]
    ZR50_L = [NaN]
    ZRmax_H = [NaN]
    ZRmax_L = [NaN]

    RfH_Zs, RfL_Zs = root_fraction_general(
        Zs, CASE_ROOT, ZR95_H, ZR50_H, ZR95_L, ZR50_L, ZRmax_H, ZRmax_L
    )

    expected_RfH_Zs = [[1.0]; zeros(15)]

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
end
