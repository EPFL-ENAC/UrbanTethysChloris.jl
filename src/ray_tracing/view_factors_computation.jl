"""
    view_factors_computation(
        XSv::Vector{FT}, YSv::Vector{FT}, dmax::FT, sz::FT, dthe::Vector{FT},
        x2::Vector{FT}, z2::Vector{FT}, x3::Vector{FT}, z3::Vector{FT},
        x4::Vector{FT}, z4::Vector{FT}, xc::FT, yc::FT, r::FT, xc2::FT,
        x5::Vector{FT}, z5::Vector{FT}
    ) where {FT<:AbstractFloat}

Compute view factors for urban canyon surfaces including trees.

# Arguments
- `XSv, YSv`: Coordinates of emitting points on emitting surface
- `x2,z2`: Coordinates of ground surface
- `x3,z3`: Coordinates of first wall surface
- `x4,z4`: Coordinates of second wall surface
- `xc,yc`: Coordinates of centre tree 1
- `xc2`: x-coordinate of centre tree 2 (yc2 = yc)
- `r`: Radius of tree
- `x5,z5`: Coordinates of sky boundary points
- `dmax`: Maximum search distance
- `sz`: Search step size
- `dthe`: Angular resolution in radians

# Returns
Tuple of view factors (VG, VW1, VW2, VS, VT1, VT2) for:
- Ground
- Wall-1
- Wall-2
- Sky
- Tree-1
- Tree-2
"""
function view_factors_computation(
    XSv::Vector{FT},
    YSv::Vector{FT},
    dmax::FT,
    sz::FT,
    dthe::Matrix{FT},
    x2::Vector{FT},
    z2::Vector{FT},
    x3::Vector{FT},
    z3::Vector{FT},
    x4::Vector{FT},
    z4::Vector{FT},
    xc::FT,
    yc::FT,
    r::FT,
    xc2::FT,
    x5::Vector{FT},
    z5::Vector{FT},
) where {FT<:AbstractFloat}
    # Initialize parameters
    spass = sqrt(FT(2)) * sz  # Pass of search
    SD = range(spass, dmax; step=spass)  # Search distance

    np = length(XSv)
    VGv = zeros(FT, np)
    VW1v = zeros(FT, np)
    VW2v = zeros(FT, np)
    VSv = zeros(FT, np)
    VT1v = zeros(FT, np)
    VT2v = zeros(FT, np)

    # For each emitting point
    for ii in 1:np
        Z = dthe[ii, :]  # Search angle [radians]
        XS = XSv[ii]
        YS = YSv[ii]
        VG = VW1 = VW2 = VS = VT1 = VT2 = zero(FT)

        # For each ray from emitting point
        for k in eachindex(Z)
            # Convert polar to cartesian coordinates for ray path
            xp, yp = pol2cart(Z[k] * ones(FT, length(SD)), collect(SD))

            # Check intersections with surfaces
            # Ground
            l1 = @SVector [x2[1], z2[1], x2[2], z2[2]]
            l2 = @SVector [XS, YS, XS + xp[end], YS + yp[end]]
            out = RayTracing.line_segment_intersect(l1, l2)
            if out.intAdjacencyMatrix
                D2 = hypot(out.intMatrixX - XS, out.intMatrixY - YS)
            else
                D2 = FT(NaN)
            end

            # Wall-1
            l1 = @SVector [x3[1], z3[1], x3[2], z3[2]]
            out = RayTracing.line_segment_intersect(l1, l2)
            if out.intAdjacencyMatrix
                D3 = hypot(out.intMatrixX - XS, out.intMatrixY - YS)
            else
                D3 = FT(NaN)
            end

            # Wall-2
            l1 = @SVector [x4[1], z4[1], x4[2], z4[2]]
            out = RayTracing.line_segment_intersect(l1, l2)
            if out.intAdjacencyMatrix
                D4 = hypot(out.intMatrixX - XS, out.intMatrixY - YS)
            else
                D4 = FT(NaN)
            end

            # Sky
            l1 = @SVector [x5[1], z5[1], x5[2], z5[2]]
            out = RayTracing.line_segment_intersect(l1, l2)
            if out.intAdjacencyMatrix
                D5 = hypot(out.intMatrixX - XS, out.intMatrixY - YS)
            else
                D5 = FT(NaN)
            end

            # Tree 1
            IC = @. (XS + xp - xc)^2 + (YS + yp - yc)^2 <= r^2
            if sum(IC) > 1
                sdi = findfirst(IC)
                DT1 = hypot(XS + xp[sdi] - XS, YS + yp[sdi] - YS)
            else
                DT1 = FT(NaN)
            end

            # Tree 2
            IC = @. (XS + xp - xc2)^2 + (YS + yp - yc)^2 <= r^2
            if sum(IC) > 1
                sdi = findfirst(IC)
                DT2 = hypot(XS + xp[sdi] - XS, YS + yp[sdi] - YS)
            else
                DT2 = FT(NaN)
            end

            # Find closest intersection
            dv = [D2, D3, D4, DT1, DT2, D5]
            dv[isnan.(dv)] .= Inf  # Set NaN to Inf
            md, pmin = findmin(dv)
            if isfinite(md)
                if pmin == 1
                    VG += 1
                elseif pmin == 2
                    VW1 += 1
                elseif pmin == 3
                    VW2 += 1
                elseif pmin == 4
                    VT1 += 1
                elseif pmin == 5
                    VT2 += 1
                elseif pmin == 6
                    VS += 1
                end
            end
        end

        # Calculate view factors for this point
        nrays = length(Z)
        VG = VG / nrays
        VW1 = VW1 / nrays
        VW2 = VW2 / nrays
        VS = VS / nrays
        VT1 = VT1 / nrays
        VT2 = VT2 / nrays

        # Normalize by total
        Sum_view = VG + VW1 + VW2 + VS + VT1 + VT2
        VGv[ii] = VG / Sum_view
        VW1v[ii] = VW1 / Sum_view
        VW2v[ii] = VW2 / Sum_view
        VSv[ii] = VS / Sum_view
        VT1v[ii] = VT1 / Sum_view
        VT2v[ii] = VT2 / Sum_view
    end

    # Calculate mean view factors
    VG = mean(VGv)
    VW1 = mean(VW1v)
    VW2 = mean(VW2v)
    VS = mean(VSv)
    VT1 = mean(VT1v)
    VT2 = mean(VT2v)

    return VG, VW1, VW2, VS, VT1, VT2
end

"""
    pol2cart(theta, r)

Convert polar coordinates to Cartesian coordinates.

# Arguments
- `theta`: Angle(s) in radians
- `r`: Radii

# Returns
Tuple of (x, y) coordinates where:
- `x = r * cos(theta)`
- `y = r * sin(theta)`
"""
function pol2cart(theta, r)
    x = r .* cos.(theta)
    y = r .* sin.(theta)
    return (x, y)
end
