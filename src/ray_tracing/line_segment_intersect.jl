"""
    line_segment_intersect(XY1::Matrix{<:Real}, XY2::Matrix{<:Real})
    line_segment_intersect(XY1::SVector{4, FT}, SVector{4, FT})

Find intersections between two sets of line segments in 2D space.

# Arguments
- `XY1`: N1×4 matrix where each row is [x1 y1 x2 y2] representing a line segment
- `XY2`: N2×4 matrix where each row is [x1 y1 x2 y2] representing a line segment

# Returns
A NamedTuple with fields:
- `intAdjacencyMatrix`: N1×N2 boolean matrix indicating intersecting segments
- `intMatrixX`: X coordinates of intersection points
- `intMatrixY`: Y coordinates of intersection points
- `intNormalizedDistance1To2`: Normalized distances from start of XY1 segments to intersections
- `intNormalizedDistance2To1`: Normalized distances from start of XY2 segments to intersections
- `parAdjacencyMatrix`: Boolean matrix indicating parallel segments
- `coincAdjacencyMatrix`: Boolean matrix indicating coincident segments
"""
function line_segment_intersect(
    XY1::AbstractMatrix{FT}, XY2::AbstractMatrix{FT}
) where {FT<:AbstractFloat}
    n_rows_1, n_cols_1 = size(XY1)
    n_rows_2, n_cols_2 = size(XY2)

    if n_cols_1 != 4 || n_cols_2 != 4
        throw(ArgumentError("Arguments must be Nx4 matrices"))
    end

    X1 = repeat(XY1[:, 1], 1, n_rows_2)
    X2 = repeat(XY1[:, 3], 1, n_rows_2)
    Y1 = repeat(XY1[:, 2], 1, n_rows_2)
    Y2 = repeat(XY1[:, 4], 1, n_rows_2)

    X3 = repeat(XY2[:, 1]', n_rows_1, 1)
    X4 = repeat(XY2[:, 3]', n_rows_1, 1)
    Y3 = repeat(XY2[:, 2]', n_rows_1, 1)
    Y4 = repeat(XY2[:, 4]', n_rows_1, 1)

    X4_X3 = X4 .- X3
    Y1_Y3 = Y1 .- Y3
    Y4_Y3 = Y4 .- Y3
    X1_X3 = X1 .- X3
    X2_X1 = X2 .- X1
    Y2_Y1 = Y2 .- Y1

    numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3
    numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3
    denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1

    u_a = numerator_a ./ denominator
    u_b = numerator_b ./ denominator

    INT_X = X1 .+ X2_X1 .* u_a
    INT_Y = Y1 .+ Y2_Y1 .* u_a
    INT_B = (u_a .>= 0) .& (u_a .<= 1) .& (u_b .>= 0) .& (u_b .<= 1)
    PAR_B = denominator .== 0
    COINC_B = (numerator_a .== 0) .& (numerator_b .== 0) .& PAR_B

    return (
        intAdjacencyMatrix=INT_B,
        intMatrixX=INT_X .* INT_B,
        intMatrixY=INT_Y .* INT_B,
        intNormalizedDistance1To2=u_a,
        intNormalizedDistance2To1=u_b,
        parAdjacencyMatrix=PAR_B,
        coincAdjacencyMatrix=COINC_B,
    )
end

function line_segment_intersect(
    XY1::SVector{4,FT}, XY2::SVector{4,FT}
) where {FT<:AbstractFloat}
    X4_X3 = XY2[3] .- XY2[1]
    Y1_Y3 = XY1[2] .- XY2[2]
    Y4_Y3 = XY2[4] .- XY2[2]
    X1_X3 = XY1[1] .- XY2[1]
    X2_X1 = XY1[3] .- XY1[1]
    Y2_Y1 = XY1[4] .- XY1[2]

    numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3
    numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3
    denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1

    u_a = numerator_a ./ denominator
    u_b = numerator_b ./ denominator

    INT_X = XY1[1] .+ X2_X1 .* u_a
    INT_Y = XY1[2] .+ Y2_Y1 .* u_a
    INT_B = (u_a .>= 0) .& (u_a .<= 1) .& (u_b .>= 0) .& (u_b .<= 1)
    PAR_B = denominator .== 0
    COINC_B = (numerator_a .== 0) .& (numerator_b .== 0) .& PAR_B

    return (
        intAdjacencyMatrix=INT_B,
        intMatrixX=INT_X .* INT_B,
        intMatrixY=INT_Y .* INT_B,
        intNormalizedDistance1To2=u_a,
        intNormalizedDistance2To1=u_b,
        parAdjacencyMatrix=PAR_B,
        coincAdjacencyMatrix=COINC_B,
    )
end
