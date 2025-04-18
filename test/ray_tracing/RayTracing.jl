using SafeTestsets

@safetestset "line_segment_intersect" begin
    include("line_segment_intersect.jl")
end

@safetestset "view_factors_computation" begin
    include("view_factors_computation.jl")
end
