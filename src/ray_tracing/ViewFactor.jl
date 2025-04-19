abstract type AbstractViewFactor{FT<:AbstractFloat} end

"""
    ViewFactorNoTrees{FT<:AbstractFloat} <: AbstractViewFactor{FT}

View factors for radiation exchange between urban surfaces without trees.

# Fields
- `F_gs`: Ground to sky view factor
- `F_gw`: Ground to wall view factor
- `F_ww`: Wall to wall view factor
- `F_wg`: Wall to ground view factor
- `F_ws`: Wall to sky view factor
- `F_sg`: Sky to ground view factor
- `F_sw`: Sky to wall view factor
"""
Base.@kwdef struct ViewFactorNoTrees{FT<:AbstractFloat} <: AbstractViewFactor{FT}
    F_gs::FT = FT(0)
    F_gw::FT = FT(0)
    F_ww::FT = FT(0)
    F_wg::FT = FT(0)
    F_ws::FT = FT(0)
    F_sg::FT = FT(0)
    F_sw::FT = FT(0)
end

"""
    ViewFactorWithTrees{FT<:AbstractFloat} <: AbstractViewFactor{FT}

View factors for radiation exchange between urban surfaces with trees.

# Fields
- `F_gs`: Ground to sky view factor
- `F_gt`: Ground to tree view factor
- `F_gw`: Ground to wall view factor
- `F_ww`: Wall to wall view factor
- `F_wt`: Wall to tree view factor
- `F_wg`: Wall to ground view factor
- `F_ws`: Wall to sky view factor
- `F_sg`: Sky to ground view factor
- `F_sw`: Sky to wall view factor
- `F_st`: Sky to tree view factor
- `F_tg`: Tree to ground view factor
- `F_tw`: Tree to wall view factor
- `F_ts`: Tree to sky view factor
- `F_tt`: Tree to tree view factor
"""
Base.@kwdef struct ViewFactorWithTrees{FT<:AbstractFloat} <: AbstractViewFactor{FT}
    F_gs::FT = FT(0)
    F_gt::FT = FT(0)
    F_gw::FT = FT(0)
    F_ww::FT = FT(0)
    F_wt::FT = FT(0)
    F_wg::FT = FT(0)
    F_ws::FT = FT(0)
    F_sg::FT = FT(0)
    F_sw::FT = FT(0)
    F_st::FT = FT(0)
    F_tg::FT = FT(0)
    F_tw::FT = FT(0)
    F_ts::FT = FT(0)
    F_tt::FT = FT(0)
end

"""
    ViewFactor{FT<:AbstractFloat} <: AbstractViewFactor{FT}

View factors for radiation exchange between urban surfaces, combining cases with and without trees.

# Fields
## Without trees (nT)
- `F_gs_nT`: Ground to sky view factor
- `F_gw_nT`: Ground to wall view factor
- `F_ww_nT`: Wall to wall view factor
- `F_wg_nT`: Wall to ground view factor
- `F_ws_nT`: Wall to sky view factor
- `F_sg_nT`: Sky to ground view factor
- `F_sw_nT`: Sky to wall view factor

## With trees (T)
- `F_gs_T`: Ground to sky view factor
- `F_gt_T`: Ground to tree view factor
- `F_gw_T`: Ground to wall view factor
- `F_ww_T`: Wall to wall view factor
- `F_wt_T`: Wall to tree view factor
- `F_wg_T`: Wall to ground view factor
- `F_ws_T`: Wall to sky view factor
- `F_sg_T`: Sky to ground view factor
- `F_sw_T`: Sky to wall view factor
- `F_st_T`: Sky to tree view factor
- `F_tg_T`: Tree to ground view factor
- `F_tw_T`: Tree to wall view factor
- `F_ts_T`: Tree to sky view factor
- `F_tt_T`: Tree to tree view factor
"""
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

function ViewFactor(
    no_trees::ViewFactorNoTrees{FT}, with_trees::ViewFactorWithTrees{FT}
) where {FT<:AbstractFloat}
    ViewFactor{FT}(
        no_trees.F_gs,
        no_trees.F_gw,
        no_trees.F_ww,
        no_trees.F_wg,
        no_trees.F_ws,
        no_trees.F_sg,
        no_trees.F_sw,
        with_trees.F_gs,
        with_trees.F_gt,
        with_trees.F_gw,
        with_trees.F_ww,
        with_trees.F_wt,
        with_trees.F_wg,
        with_trees.F_ws,
        with_trees.F_sg,
        with_trees.F_sw,
        with_trees.F_st,
        with_trees.F_tg,
        with_trees.F_tw,
        with_trees.F_ts,
        with_trees.F_tt,
    )
end

"""
    ViewFactorPoint{FT<:AbstractFloat} <: AbstractViewFactor{FT}

View factors from a single point to urban surfaces.

# Fields
- `F_pg`: Point to ground view factor
- `F_ps`: Point to sky view factor
- `F_pt`: Point to tree view factor
- `F_pwLeft`: Point to left wall view factor
- `F_pwRight`: Point to right wall view factor
"""
Base.@kwdef struct ViewFactorPoint{FT<:AbstractFloat} <: AbstractViewFactor{FT}
    F_pg::FT
    F_ps::FT
    F_pt::FT
    F_pwLeft::FT
    F_pwRight::FT
end
