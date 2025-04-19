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

    vf = ViewFactor{FT}(;
        F_gs_T=view_factor[3, 6],
        F_gt_T=view_factor[3, 4] + view_factor[3, 5],
        F_gw_T=(view_factor[3, 1] + view_factor[3, 2])/2,
        F_ww_T=view_factor[1, 2],
        F_wt_T=view_factor[1, 4] + view_factor[1, 5],
        F_wg_T=view_factor[1, 3],
        F_ws_T=view_factor[1, 6],
        F_ts_T=view_factor[4, 6] * a_mask,
        F_tw_T=(view_factor[4, 1] + view_factor[4, 2])/2 * a_mask,
        F_tt_T=view_factor[4, 5] * a_mask,
        F_tg_T=view_factor[4, 3] * a_mask,
        F_sg_T=view_factor[6, 3],
        F_sw_T=(view_factor[6, 1] + view_factor[6, 2])/2,
        F_st_T=view_factor[6, 4] + view_factor[6, 5],
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
