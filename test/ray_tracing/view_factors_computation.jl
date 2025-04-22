using Test
using UrbanTethysChloris.RayTracing: view_factors_computation

FT = Float64

@testset "view_factors_computation" begin
    dmax = 1.517000264932407
    dthe = repeat(LinRange(-pi / 2, pi / 2, 10)', 10, 1)
    r = 0.05
    sz = 1e-3
    x2 = FT[1, 2]
    x3 = FT[1, 1]
    x4 = FT[2, 2]
    x5 = FT[1, 2]
    xc = 1.222413793103448
    xc2 = 1.777586206896552
    XSv = fill(1.0000000001, 10)
    yc = 0.812068965517241
    YSv = [
        0.972604932176721,
        0.605318590650978,
        0.338236165549055,
        0.927983722291575,
        0.898424694301282,
        0.850706460713442,
        0.256791917068097,
        0.285495500891502,
        0.779946910428948,
        0.701395228244418,
    ]
    z2 = FT[0, 0]
    z3 = [1.120689655172414, 0]
    z4 = [0, 1.120689655172414]
    z5 = z3 + z4

    @testset "MATLAB - Matrix" begin
        VG, VW1, VW2, VS, VT1, VT2 = view_factors_computation(
            XSv, YSv, dmax, sz, dthe, x2, z2, x3, z3, x4, z4, xc, yc, r, xc2, x5, z5
        )

        @test VG ≈ 0.33
        @test VW1 == 0
        @test VW2 ≈ 0.19
        @test VS ≈ 0.33
        @test VT1 ≈ 0.11
        @test VT2 ≈ 0.04
    end

    @testset "MATLAB - Vector" begin
        VG, VW1, VW2, VS, VT1, VT2 = view_factors_computation(
            XSv, YSv, dmax, sz, dthe[1, :], x2, z2, x3, z3, x4, z4, xc, yc, r, xc2, x5, z5
        )

        @test VG ≈ 0.33
        @test VW1 == 0
        @test VW2 ≈ 0.19
        @test VS ≈ 0.33
        @test VT1 ≈ 0.11
        @test VT2 ≈ 0.04
    end
end
