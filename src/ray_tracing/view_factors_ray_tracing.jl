"""
    view_factors_ray_tracing(
        H::FT,
        W::FT,
        a::FT,
        ht::FT,
        d::FT,
        person::ModelComponents.Parameters.PersonParameters{FT},
        mc_sample_size::Int,
        n_rays::Int,
    ) where {FT<:AbstractFloat}

Calculate view factors for urban canyon surfaces using ray tracing method.

# Arguments
- `H`: canyon height [m]
- `W`: canyon width [m]
- `a`: normalized tree radius [-]
- `ht`: normalized tree height [-]
- `d`: normalized tree distance from wall [-]
- `person`: PersonParameters object containing position information
- `mc_sample_size`: number of Monte Carlo samples for ray tracing
- `n_rays`: number of rays to emit per sample point

# Returns
Tuple containing:
- `vf`: ViewFactorWithTrees object containing view factors between surfaces
- `vfp`: ViewFactorPoint object containing point-specific view factors
"""
function view_factors_ray_tracing(
    H::FT,
    W::FT,
    a::FT,
    ht::FT,
    d::FT,
    person::ModelComponents.Parameters.PersonParameters{FT},
    mc_sample_size::Int,
    n_rays::Int,
) where {FT<:AbstractFloat}
    view_factor = zeros(FT, 7, 6)

    # Compute view factors for each surface
    for option_surface in 1:7
        VG, VW1, VW2, VS, VT1, VT2 = view_factors_geometry(
            H, W, a, ht, d, person, option_surface, mc_sample_size, n_rays
        )
        view_factor[option_surface, 1] = VW1  # towards wall 1
        view_factor[option_surface, 2] = VW2  # towards wall 2
        view_factor[option_surface, 3] = VG   # towards ground
        view_factor[option_surface, 4] = VT1  # towards tree 1
        view_factor[option_surface, 5] = VT2  # towards tree 2
        view_factor[option_surface, 6] = VS   # towards sky
    end

    # Check view factors sum to 1
    svf_test = sum(view_factor; dims=2)
    any(x -> !(0.9999 ≤ x ≤ 1.0001), svf_test) &&
        @warn "View factors do not add up to 1 after ray tracing"

    # Eliminate self-view factors and rescale
    for i in 1:6
        view_factor[i, :] ./= (1 - view_factor[i, i])
        view_factor[i, i] = 0
    end

    # View factor assignments
    a_mask = a > 0

    vf = ViewFactorWithTrees{FT}(;
        F_gs=view_factor[3, 6],
        F_gt=view_factor[3, 4] + view_factor[3, 5],
        F_gw=(view_factor[3, 1] + view_factor[3, 2])/2,
        F_ww=view_factor[1, 2],
        F_wt=view_factor[1, 4] + view_factor[1, 5],
        F_wg=view_factor[1, 3],
        F_ws=view_factor[1, 6],
        F_ts=view_factor[4, 6] * a_mask,
        F_tw=(view_factor[4, 1] + view_factor[4, 2])/2 * a_mask,
        F_tt=view_factor[4, 5] * a_mask,
        F_tg=view_factor[4, 3] * a_mask,
        F_sg=view_factor[6, 3],
        F_sw=(view_factor[6, 1] + view_factor[6, 2])/2,
        F_st=view_factor[6, 4] + view_factor[6, 5],
    )

    vfp = ViewFactorPoint{FT}(;
        F_pg=view_factor[7, 3],
        F_ps=view_factor[7, 6],
        F_pt=view_factor[7, 4] + view_factor[7, 5],
        F_pwLeft=view_factor[7, 1],
        F_pwRight=view_factor[7, 2],
    )

    return vf, vfp
end

"""
    view_factors_ray_tracing_reciprocity(H::FT, W::FT, a::FT, ht::FT, d::FT, person::PersonParameters{FT},
                           mc_sample_size::Int, n_rays::Int) where {FT<:AbstractFloat}

Compute reciprocal view factors for urban canyon geometry.

# Arguments
- `H`: canyon height [m]
- `W`: canyon width [m]
- `a`: normalized tree radius [-]
- `ht`: normalized tree height [-]
- `d`: normalized tree distance from wall [-]
- `person`: PersonParameters object containing position information
- `mc_sample_size`: number of emitting points per surface
- `n_rays`: number of rays emitted per point

# Returns
Tuple containing (ViewFactorWithTrees, ViewFactorPoint, ViewFactorWithTrees)
"""
function view_factors_ray_tracing_reciprocity(
    H::FT,
    W::FT,
    a::FT,
    ht::FT,
    d::FT,
    person::ModelComponents.Parameters.PersonParameters{FT},
    mc_sample_size::Int,
    n_rays::Int,
) where {FT<:AbstractFloat}

    # Get raw view factors from ray tracing
    vf_raw, vf_point = view_factors_ray_tracing(
        H, W, a, ht, d, person, mc_sample_size, n_rays
    )

    # Normalize dimensions
    h = H/W
    w = W/W

    # Point view factors already correctly populated by view_factors_ray_tracing

    # Compute reciprocal view factors
    if a == 0
        # Case without trees
        F_gs_T = vf_raw.F_gs
        F_gw_T = (1 - F_gs_T)/2  # factor 1/2 because there are 2 walls seen by ground
        F_gt_T = zero(FT)

        F_sg_T = F_gs_T * w/w
        F_sw_T = F_gw_T * w/w
        F_st_T = zero(FT)

        F_wg_T = F_gw_T * w/h
        F_ws_T = F_sw_T * w/h
        F_ww_T = 1 - F_wg_T - F_ws_T
        F_wt_T = zero(FT)

        F_tg_T = F_ts_T = F_tw_T = F_tt_T = zero(FT)

        sum = zeros(FT, 4)
        sum[1] = F_gs_T + 2*F_gw_T
        sum[2] = F_ww_T + F_wg_T + F_ws_T
        sum[3] = F_sg_T + 2*F_sw_T
        sum[4] = zero(FT)

        sum2 = zeros(FT, 4)
        sum2[1] = F_sg_T*w/w + 2*F_wg_T*h/w
        sum2[2] = F_ww_T*h/h + F_gw_T*w/h + F_sw_T*w/h
        sum2[3] = F_gs_T*w/w + 2*F_ws_T*h/w
        sum2[4] = zero(FT)
    else
        # Case with trees
        Atree = 2 * 2π * a

        F_gs_T = vf_raw.F_gs
        F_gt_T = vf_raw.F_gt
        F_gw_T = (1 - F_gs_T - F_gt_T)/2

        F_sg_T = F_gs_T * w/w
        F_st_T = vf_raw.F_st
        F_sw_T = (1 - F_sg_T - F_st_T)/2

        F_wg_T = F_gw_T * w/h
        F_ws_T = F_sw_T * w/h
        F_wt_T = vf_raw.F_wt
        F_ww_T = 1 - F_wg_T - F_ws_T - F_wt_T

        F_ts_T = F_st_T * w/Atree
        F_tw_T = F_wt_T * h/Atree
        F_tg_T = F_gt_T * w/Atree
        F_tt_T = 1 - F_ts_T - 2*F_tw_T - F_tg_T

        sum = zeros(FT, 4)
        sum[1] = F_gs_T + 2*F_gw_T + F_gt_T
        sum[2] = F_ww_T + F_wg_T + F_ws_T + F_wt_T
        sum[3] = F_sg_T + 2*F_sw_T + F_st_T
        sum[4] = F_tg_T + 2*F_tw_T + F_ts_T + F_tt_T

        sum2 = zeros(FT, 4)
        sum2[1] = F_sg_T*w/w + 2*F_wg_T*h/w + F_tg_T*Atree/w
        sum2[2] = F_ww_T*h/h + F_gw_T*w/h + F_sw_T*w/h + F_tw_T*Atree/h
        sum2[3] = F_gs_T*w/w + 2*F_ws_T*h/w + F_ts_T*Atree/w
        sum2[4] = F_gt_T*w/Atree + 2*F_wt_T*h/Atree + F_st_T*w/Atree + F_tt_T*Atree/Atree
    end

    # Check reciprocity conditions
    if a > 0
        for sum_i in sum
            if !(0.9999 ≤ sum_i ≤ 1.0001)
                @warn "View factors do not sum to 1 (sum = $sum_i)"
            end
        end
    else
        for sum_i in view(sum, 1:3)
            if !(0.9999 ≤ sum_i ≤ 1.0001)
                @warn "View factors do not sum to 1 (sum = $sum_i)"
            end
        end
    end

    # Create view factor structs
    vf = ViewFactorWithTrees{FT}(;
        F_gs=F_gs_T,
        F_gt=F_gt_T,
        F_gw=F_gw_T,
        F_ww=F_ww_T,
        F_wt=F_wt_T,
        F_wg=F_wg_T,
        F_ws=F_ws_T,
        F_sg=F_sg_T,
        F_sw=F_sw_T,
        F_st=F_st_T,
        F_tg=F_tg_T,
        F_tw=F_tw_T,
        F_ts=F_ts_T,
        F_tt=F_tt_T,
    )

    return vf, vf_point, vf_raw
end
