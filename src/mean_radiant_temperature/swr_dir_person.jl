"""
    swr_dir_person(
        SWR_dir::FT,
        zeta_S::FT,
        theta_Z::FT,
        BooleanInSun::FT,
    ) where {FT<:AbstractFloat}

Calculate direct shortwave radiation onto person.

# Arguments
- `SWR_dir`: Direct shortwave radiation [W/m²]
- `zeta_S`: Solar azimuth angle [rad]
- `theta_Z`: Solar zenith angle [rad]
- `BooleanInSun`: Person is in sun, partially, or in shade (0=shade, 1=sun, 0-1=partial)

# Returns
- `SWRdir_Person::FT`: Direct shortwave radiation on person [W/m²]
- `SWRdir_in_top::FT`: Direct radiation on top surface [W/m²]
- `SWRdir_in_bottom::FT`: Direct radiation on bottom surface [W/m²]
- `SWRdir_in_east::FT`: Direct radiation on east surface [W/m²]
- `SWRdir_in_south::FT`: Direct radiation on south surface [W/m²]
- `SWRdir_in_west::FT`: Direct radiation on west surface [W/m²]
- `SWRdir_in_north::FT`: Direct radiation on north surface [W/m²]
"""
function swr_dir_person(
    SWR_dir::FT, zeta_S::FT, theta_Z::FT, BooleanInSun::FT
) where {FT<:AbstractFloat}
    if zeta_S < 0
        @warn "Solar azimuth angle smaller than 0 which is not possible in this formulation. Please change to 0 to 2pi"
    end

    # Top and bottom direct solar radiation
    SWRdir_in_top = SWR_dir * BooleanInSun
    SWRdir_in_bottom = zero(FT)

    # Solar radiation on east component
    if 0 < zeta_S < π
        SWRdir_in_east = SWR_dir * BooleanInSun * tan(theta_Z) * abs(sin(zeta_S))
    else
        SWRdir_in_east = zero(FT)
    end

    # Solar radiation on south component
    if π/2 < zeta_S < 3π/2
        SWRdir_in_south = SWR_dir * BooleanInSun * tan(theta_Z) * abs(sin(zeta_S - π/2))
    else
        SWRdir_in_south = zero(FT)
    end

    # Solar radiation on west component
    if π < zeta_S < 2π
        SWRdir_in_west = SWR_dir * BooleanInSun * tan(theta_Z) * abs(sin(zeta_S - π))
    else
        SWRdir_in_west = zero(FT)
    end

    # Solar radiation on north component
    if zeta_S > 3π/2 || zeta_S < π/2
        SWRdir_in_north = SWR_dir * BooleanInSun * tan(theta_Z) * abs(sin(zeta_S - 3π/2))
    else
        SWRdir_in_north = zero(FT)
    end

    # Total direct solar radiation arriving onto human
    SWRdir_Person =
        FT(0.06) * (SWRdir_in_top + SWRdir_in_bottom) +
        FT(0.22) * (SWRdir_in_north + SWRdir_in_east + SWRdir_in_south + SWRdir_in_west)

    return SWRdir_Person,
    SWRdir_in_top, SWRdir_in_bottom, SWRdir_in_east, SWRdir_in_south, SWRdir_in_west,
    SWRdir_in_north
end
