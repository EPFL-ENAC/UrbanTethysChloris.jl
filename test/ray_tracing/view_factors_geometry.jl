using Test
using UrbanTethysChloris.RayTracing: view_factors_geometry
using UrbanTethysChloris.ModelComponents.Parameters: PersonParameters
using Random

FT = Float64

@testset "view_factors_geometry" begin

    # Test inputs
    a = FT(0.05)
    d = FT(0.222413793103448)
    H = FT(6.5)
    ht = FT(0.812068965517241)
    mc_sample_size = 10
    n_rays = 10
    W = FT(5.8)

    person = PersonParameters{FT}(
        PositionPx=W/2, PositionPz=1.1, PersonWidth=0.03, PersonHeight=0.11, HeightWind=1.1
    )

    # Test each option_surface value
    for option_surface in 1:7
        VG, VW1, VW2, VS, VT1, VT2 = view_factors_geometry(
            H,
            W,
            option_surface == 1 ? zero(FT) : a,
            ht,
            d,
            person,
            option_surface,
            mc_sample_size,
            n_rays,
        )

        # Basic sanity checks
        @test 0 ≤ VG ≤ 1
        @test 0 ≤ VW1 ≤ 1
        @test 0 ≤ VW2 ≤ 1
        @test 0 ≤ VS ≤ 1
        @test 0 ≤ VT1 ≤ 1
        @test 0 ≤ VT2 ≤ 1

        # Sum should be approximately 1
        @test isapprox(VG + VW1 + VW2 + VS + VT1 + VT2, 1.0, atol=1e-10)
    end
end

@testset "No trees" begin

    # Test inputs
    a = 0.0
    d = 0.0
    H = FT(6.5)
    ht = -0.00065
    mc_sample_size = 10
    n_rays = 10
    W = FT(5.8)

    person = PersonParameters{FT}(
        PositionPx=W/2, PositionPz=1.1, PersonWidth=0.03, PersonHeight=0.11, HeightWind=1.1
    )

    @testset "NaNs" begin
        # Test each option_surface value
        Random.seed!(1234)  # Set a random seed for reproducibility
        VG, VW1, VW2, VS, VT1, VT2 = view_factors_geometry(
            H, W, a, ht, d, person, 4, mc_sample_size, n_rays
        )

        # Basic sanity checks
        @test isnan(VG)
        @test isnan(VW1)
        @test isnan(VW2)
        @test isnan(VS)
        @test isnan(VT1)
        @test isnan(VT2)
    end

    @testset "No NaNs" begin
        Random.seed!(123)  # Set a random seed for reproducibility
        VG, VW1, VW2, VS, VT1, VT2 = view_factors_geometry(
            H, W, a, ht, d, person, 5, mc_sample_size, n_rays
        )

        # Basic sanity checks
        @test 0 ≤ VG ≤ 1
        @test 0 ≤ VW1 ≤ 1
        @test 0 ≤ VW2 ≤ 1
        @test 0 ≤ VS ≤ 1
        @test 0 ≤ VT1 ≤ 1
        @test 0 ≤ VT2 ≤ 1

        # Sum should be approximately 1
        @test isapprox(VG + VW1 + VW2 + VS + VT1 + VT2, 1.0, atol=1e-10)
    end
end
