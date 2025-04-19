using SafeTestsets

@safetestset "line_segment_intersect" begin
    include("line_segment_intersect.jl")
end

@safetestset "view_factors_computation" begin
    include("view_factors_computation.jl")
end

@safetestset "view_factors_geometry" begin
    include("view_factors_geometry.jl")
end

@safetestset "view_factors_ray_tracing" begin
    include("view_factors_ray_tracing.jl")
end

@safetestset "view_factors_analytical" begin
    include("view_factors_analytical.jl")
end
