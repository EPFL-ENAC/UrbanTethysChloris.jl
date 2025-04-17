using Test
using UrbanTethysChloris.RayTracing: line_segment_intersect

FT = Float64

@testset "MATLAB" begin
    XY1 = FT[0.1 0 2.3 0; 1.5 0 2.7 1.8]
    XY2 = FT[0 0 1.3 1.4; 0 1.6 1.7 0]

    out = line_segment_intersect(XY1, XY2)

    @test out.intAdjacencyMatrix == [false true; false true]
    @test out.intMatrixX ≈ [0 1.7; 0 1.577108433734940]
    @test out.intMatrixY ≈ [0 0; 0 0.115662650602410]
    @test out.intNormalizedDistance1To2 ≈ [
        -0.045454545454545 0.727272727272727;
        3.181818181818181 0.064257028112450
    ]
    @test out.intNormalizedDistance2To1 ≈ [
        0 1;
        4.090909090909090 0.927710843373494
    ]
    @test out.parAdjacencyMatrix == [false false; false false]
    @test out.coincAdjacencyMatrix == [false false; false false]
end

@testset "Wrong dimensions" begin
    XY1 = FT[0.1 0 2.3 0; 1.5 0 2.7 1.8]
    XY2 = FT[0 1.3 1.4; 0 1.6 1.7]

    @test_throws ArgumentError line_segment_intersect(XY1, XY2[:, 1:3])
    @test_throws ArgumentError line_segment_intersect(XY1[:, 1:3], XY2)
end
