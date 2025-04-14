abstract type AbstractViewFactor{FT<:AbstractFloat} end

Base.@kwdef struct ViewFactor{FT<:AbstractFloat} <: AbstractViewFactor{FT}
    F_gs_nT::FT = FT(0)
    F_gw_nT::FT = FT(0)
    F_ww_nT::FT = FT(0)
    F_wg_nT::FT = FT(0)
    F_ws_nT::FT = FT(0)
    F_sg_nT::FT = FT(0)
    F_sw_nT::FT = FT(0)
    F_gs_T::FT = FT(0)
    F_gt_T::FT = FT(0)
    F_gw_T::FT = FT(0)
    F_ww_T::FT = FT(0)
    F_wt_T::FT = FT(0)
    F_wg_T::FT = FT(0)
    F_ws_T::FT = FT(0)
    F_sg_T::FT = FT(0)
    F_sw_T::FT = FT(0)
    F_st_T::FT = FT(0)
    F_tg_T::FT = FT(0)
    F_tw_T::FT = FT(0)
    F_ts_T::FT = FT(0)
    F_tt_T::FT = FT(0)
end

Base.@kwdef struct ViewFactorPoint{FT<:AbstractFloat} <: AbstractViewFactor{FT}
    F_pg::FT
    F_ps::FT
    F_pt::FT
    F_pwLeft::FT
    F_pwRight::FT
end
