"""
    swr_diff_person(
        SWRout_t::NamedTuple,
        LWRout_t::NamedTuple,
        MeteoData::NamedTuple,
        ViewFactorPoint::RayTracing.ViewFactorPoint{FT},
        TimeOfMaxSolAlt::FT,
        TimeHr::FT,
    ) where {FT<:AbstractFloat}

Calculate diffuse shortwave and longwave radiation onto point.

# Arguments
- `SWRout_t`: Outgoing shortwave radiation from surfaces [W/m²]
- `LWRout_t`: Outgoing longwave radiation from surfaces [W/m²]
- `MeteoData`: Meteorological data
- `ViewFactorPoint`: View factors from point to surfaces
- `TimeOfMaxSolAlt`: Time of maximum solar altitude [h]
- `TimeHr`: Current time [h]

# Returns
- `SWRdiff_Person::FT`: Diffuse shortwave radiation at point [W/m²]
- `LWR_Person::FT`: Longwave radiation at point [W/m²]
"""
function swr_diff_person(
    SWRout_t::Radiation.RadiationFluxes{FT},
    LWRout_t::Radiation.RadiationFluxes{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    ViewFactorPoint::ViewFactorPoint{FT},
    TimeOfMaxSolAlt::FT,
    TimeHr::FT,
) where {FT<:AbstractFloat}
    if TimeHr ≤ TimeOfMaxSolAlt
        SWRdiff_Person =
            ViewFactorPoint.F_pg * SWRout_t.TotalGround +
            ViewFactorPoint.F_ps * MeteoData.SW_diff +
            ViewFactorPoint.F_pt * SWRout_t.Tree +
            ViewFactorPoint.F_pwLeft * SWRout_t.WallSun +
            ViewFactorPoint.F_pwRight * SWRout_t.WallShade

        LWR_Person =
            ViewFactorPoint.F_pg * LWRout_t.TotalGround +
            ViewFactorPoint.F_ps * MeteoData.LWR +
            ViewFactorPoint.F_pt * LWRout_t.Tree +
            ViewFactorPoint.F_pwLeft * LWRout_t.WallSun +
            ViewFactorPoint.F_pwRight * LWRout_t.WallShade
    else
        SWRdiff_Person =
            ViewFactorPoint.F_pg * SWRout_t.TotalGround +
            ViewFactorPoint.F_ps * MeteoData.SW_diff +
            ViewFactorPoint.F_pt * SWRout_t.Tree +
            ViewFactorPoint.F_pwRight * SWRout_t.WallSun +
            ViewFactorPoint.F_pwLeft * SWRout_t.WallShade

        LWR_Person =
            ViewFactorPoint.F_pg * LWRout_t.TotalGround +
            ViewFactorPoint.F_ps * MeteoData.LWR +
            ViewFactorPoint.F_pt * LWRout_t.Tree +
            ViewFactorPoint.F_pwRight * LWRout_t.WallSun +
            ViewFactorPoint.F_pwLeft * LWRout_t.WallShade
    end

    return SWRdiff_Person, LWR_Person
end

function swr_diff_person(
    SWRout_t::NamedTuple,
    LWRout_t::NamedTuple,
    MeteoData::NamedTuple,
    ViewFactorPoint::ViewFactorPoint{FT},
    TimeOfMaxSolAlt::FT,
    TimeHr::FT,
) where {FT<:AbstractFloat}
    if TimeHr <= TimeOfMaxSolAlt
        SWRdiff_Person =
            ViewFactorPoint.F_pg * SWRout_t.SWRoutTotalGround +
            ViewFactorPoint.F_ps * MeteoData.SW_diff +
            ViewFactorPoint.F_pt * SWRout_t.SWRoutTree +
            ViewFactorPoint.F_pwLeft * SWRout_t.SWRoutWallSun +
            ViewFactorPoint.F_pwRight * SWRout_t.SWRoutWallShade

        LWR_Person =
            ViewFactorPoint.F_pg * LWRout_t.LWRoutTotalGround +
            ViewFactorPoint.F_ps * MeteoData.LWR +
            ViewFactorPoint.F_pt * LWRout_t.LWRoutTree +
            ViewFactorPoint.F_pwLeft * LWRout_t.LWRoutWallSun +
            ViewFactorPoint.F_pwRight * LWRout_t.LWRoutWallShade
    else
        SWRdiff_Person =
            ViewFactorPoint.F_pg * SWRout_t.SWRoutTotalGround +
            ViewFactorPoint.F_ps * MeteoData.SW_diff +
            ViewFactorPoint.F_pt * SWRout_t.SWRoutTree +
            ViewFactorPoint.F_pwRight * SWRout_t.SWRoutWallSun +
            ViewFactorPoint.F_pwLeft * SWRout_t.SWRoutWallShade

        LWR_Person =
            ViewFactorPoint.F_pg * LWRout_t.LWRoutTotalGround +
            ViewFactorPoint.F_ps * MeteoData.LWR +
            ViewFactorPoint.F_pt * LWRout_t.LWRoutTree +
            ViewFactorPoint.F_pwRight * LWRout_t.LWRoutWallSun +
            ViewFactorPoint.F_pwLeft * LWRout_t.LWRoutWallShade
    end

    return SWRdiff_Person, LWR_Person
end
