"""
    view_factors_analytical(H::FT, W::FT) where {FT<:AbstractFloat}

Compute analytical view factors without trees according to Harman et al. 2004

# Arguments
- `H`: canyon height [m]
- `W`: canyon width [m]

# Returns
- `ViewFactorNoTrees`: struct containing view factors
"""
function view_factors_analytical(H::FT, W::FT) where {FT<:AbstractFloat}
    ratio = H/W

    F_gs = sqrt(1 + ratio^2) - ratio
    F_gw = FT(0.5) * (1 - F_gs)  # factor 0.5 because there are 2 walls seen by ground

    F_ww = sqrt(1 + (1/ratio)^2) - 1/ratio
    F_wg = FT(0.5) * (1 - F_ww)
    F_ws = FT(0.5) * (1 - F_ww)

    F_sg = F_gs
    F_sw = ratio * F_ws

    return ViewFactorNoTrees{FT}(;
        F_gs=F_gs, F_gw=F_gw, F_ww=F_ww, F_wg=F_wg, F_ws=F_ws, F_sg=F_sg, F_sw=F_sw
    )
end
