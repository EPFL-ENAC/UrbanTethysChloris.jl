"""
    view_factor_internal(
        Hbuild::FT,
        Wroof::FT
    ) where {FT<:AbstractFloat}

Calculate view factors for internal building surfaces.

# Arguments
- `Hbuild`: Building/canyon height [m]
- `Wroof`: Roof/canyon width [m]

# Returns
- `ViewFactor`: View factors between internal surfaces
"""
function view_factor_internal(Hbuild::FT, Wroof::FT) where {FT<:AbstractFloat}
    ratio = Hbuild / Wroof

    # Calculate view factors
    F_gc = sqrt(1 + ratio^2) - ratio
    F_gw = (1 - F_gc) / 2  # factor 0.5 because of 4 walls (2 external, 2 internal)

    F_ww = sqrt(1 + (1/ratio)^2) - 1/ratio
    F_wg = (1 - F_ww) / 2
    F_wc = (1 - F_ww) / 2

    F_cg = F_gc
    F_cw = ratio * F_wc

    # Check unity of view factor sums
    h = Hbuild / Wroof
    w = one(FT)

    Sum_g = F_gc + F_gw * 2
    Sum_w = F_ww + F_wg + F_wc
    Sum_c = F_cg + 2 * F_cw

    Sum_g2 = F_wg * h/w * 2 + F_cg
    Sum_w2 = F_gw * w/h + F_ww + F_cw * w/h
    Sum_c2 = F_gc + 2 * F_wc * h/w

    # Return named tuple of view factors
    ViewFactor = ViewFactorInternal(;
        F_gc=F_gc, F_gw=F_gw, F_ww=F_ww, F_wg=F_wg, F_wc=F_wc, F_cg=F_cg, F_cw=F_cw
    )

    return F_gc, F_gw, F_ww, F_wg, F_wc, F_cg, F_cw, ViewFactor
end
