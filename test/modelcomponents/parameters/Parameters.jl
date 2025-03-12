using SafeTestsets

@safetestset "SurfaceFractions" begin
    include("SurfaceFractions.jl")
end

@safetestset "UrbanGeometryParameters" begin
    include("UrbanGeometryParameters.jl")
end
