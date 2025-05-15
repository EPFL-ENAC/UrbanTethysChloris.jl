using Test
using UrbanTethysChloris.RayTracing:
    view_factors_ray_tracing, view_factors_ray_tracing_reciprocity
using UrbanTethysChloris.ModelComponents.Parameters: PersonParameters
using Random

FT = Float64

a = FT(0.05)
d = FT(0.222413793103448)
H = FT(6.5)
ht = FT(0.812068965517241)
mc_sample_size = 10
n_rays = 10
W = FT(5.8)

person = PersonParameters{FT}(;
    PositionPx=W / 2, PositionPz=1.1, PersonWidth=0.03, PersonHeight=0.11, HeightWind=1.1
)

@testset "view_factors_ray_tracing" begin
    @testset "Trees" begin
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

    @testset "No trees" begin
        Random.seed!(123)  # Set seed for reproducibility
        view_factor, vfp = view_factors_ray_tracing(
            H, W, zero(FT), ht, d, person, mc_sample_size, n_rays
        )

        # Check view factors are within expected ranges
        @test all(
            0 ≤ getfield(view_factor, field) ≤ 1 for
            field in fieldnames(typeof(view_factor))
        )
        @test all(0 ≤ getfield(vfp, field) ≤ 1 for field in fieldnames(typeof(vfp)))
    end
end

@testset "view_factors_ray_tracing_reciprocity" begin
    Random.seed!(123)
    vf, vfp, vf_raw = view_factors_ray_tracing_reciprocity(
        H, W, a, ht, d, person, mc_sample_size, n_rays
    )

    # Test view factor ranges
    @test all(0 ≤ getfield(vf, field) ≤ 1 for field in fieldnames(typeof(vf)))
    @test all(0 ≤ getfield(vfp, field) ≤ 1 for field in fieldnames(typeof(vfp)))
    @test all(0 ≤ getfield(vf_raw, field) ≤ 1 for field in fieldnames(typeof(vf_raw)))
end
