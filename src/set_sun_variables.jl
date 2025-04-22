"""
    set_sun_variables(datam, deltagmt, lon, lat, t_bef, t_aft)

Calculate solar position variables.

Returns tuple of (h_S, zeta_S, T_sunrise, T_sunset) where:
- h_S: solar altitude [rad]
- zeta_S: Sun's azimuth [rad]
- T_sunrise: sunrise time [h]
- T_sunset: sunset time [h]
"""
function set_sun_variables(
    datam::Vector{Int}, deltagmt::FT, lon::FT, lat::FT, t_bef::FT, t_aft::FT
) where {FT<:AbstractFloat}

    # Calculate Julian day
    jday = dayofyear(datam[1], datam[2], datam[3])
    nowhr = datam[4] + datam[5]/60 + datam[6]/3600

    # Solar declination
    delta_s = 23.45 * π/180 * cos(2π/365 * (172 - jday))

    # Time difference between standard and local meridian
    delta_tsl = if lon < 0
        -1/15 * (15*abs(deltagmt) - abs(lon))
    else
        1/15 * (15*abs(deltagmt) - abs(lon))
    end

    # Generate time array
    t = range(nowhr-t_bef, nowhr+t_aft; step=0.0166666)

    # Calculate hour angle
    tau_s = [
        if ti < (12 + delta_tsl)
            15π/180*(ti + 12 - delta_tsl)
        else
            15π/180*(ti - 12 - delta_tsl)
        end for ti in t
    ]

    # Solar altitude
    lat_rad = lat * π/180
    sinh_s = sin(lat_rad)*sin(delta_s) .+ cos(lat_rad)*cos(delta_s)*cos.(tau_s)
    h_s = mean(asin.(sinh_s))

    # Sun's azimuth
    zeta_s = atan.(-sin.(tau_s) ./ (tan(delta_s)*cos(lat_rad) .- sin(lat_rad)*cos.(tau_s)))

    # Adjust azimuth based on hour angle
    for i in eachindex(tau_s)
        if 0 < tau_s[i] <= π
            zeta_s[i] += zeta_s[i] > 0 ? π : 2π
        elseif π <= tau_s[i] <= 2π && zeta_s[i] < 0
            zeta_s[i] += π
        end
    end
    zeta_s = mean(zeta_s)

    # Sunrise and sunset times
    t_sunrise = 180/(15π) * (2π - acos(-tan(delta_s)*tan(lat_rad))) - 12
    t_sunset = 180/(15π) * acos(-tan(delta_s)*tan(lat_rad)) + 12

    return h_s, real(zeta_s), real(t_sunrise), real(t_sunset)
end
