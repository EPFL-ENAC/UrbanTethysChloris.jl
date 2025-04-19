using Test
using UrbanTethysChloris.RayTracing: view_factors_ray_tracing
using UrbanTethysChloris.ModelComponents.Parameters: PersonParameters

FT = Float64

@testset "view_factors_ray_tracing" begin
    a = FT(0.05)
    d = FT(0.222413793103448)
    H = FT(6.5)
    ht = FT(0.812068965517241)
    mc_sample_size = 10
    n_rays = 10
    W = FT(5.8)

    person = PersonParameters{FT}(;
        PositionPx=W / 2,
        PositionPz=1.1,
        PersonWidth=0.03,
        PersonHeight=0.11,
        HeightWind=1.1,
    )

    @testset "MATLAB" begin
        view_factor, vfp = view_factors_ray_tracing(
            H, W, a, ht, d, person, mc_sample_size, n_rays
        )

        # Check view factors are within expected ranges
        @test all(
            0 ≤ getfield(view_factor, field) ≤ 1 for
            field in fieldnames(typeof(view_factor))
        )
        @test all(0 ≤ getfield(vfp, field) ≤ 1 for field in fieldnames(typeof(vfp)))
    end
end
