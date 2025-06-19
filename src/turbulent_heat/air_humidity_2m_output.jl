"""
    air_humidity_2m_output(
        q2m::FT, T2m::FT, Timp::FT, Tbare::FT, Tveg::FT, Tcan::FT, qcan::FT,
        rap_can2m::FT, rap_can2m_Inv::FT, rb_L::FT, alp_soil_bare::FT, r_soil_bare::FT,
        alp_soil_veg::FT, r_soil_veg::FT, rs_sun_L::FT, rs_shd_L::FT,
        dw_L::FT, Fsun_L::FT, Fshd_L::FT, FractionsGround, ParVegGround,
        Eimp::FT, Ebare::FT, Eveg_int::FT, Eveg_pond::FT, Eveg_soil::FT, TEveg::FT, Pre::FT,
        Humidity_ittm, fconv::FT, MeteoData, Gemeotry_m, rho_atm::FT, Zp1::FT, ParCalculation
    ) where {FT<:AbstractFloat}

Calculate the air humidity at 2m level and compute vapor flux differences with additional outputs.

# Arguments
- `q2m`: specific humidity at 2m [kg/kg]
- `T2m`: temperature at 2m [K]
- `Timp`: impervious surface temperature [K]
- `Tbare`: bare soil temperature [K]
- `Tveg`: vegetation temperature [K]
- `Tcan`: canyon air temperature [K]
- `qcan`: canyon specific humidity [kg/kg]
[...remaining arguments same as air_humidity_2m...]

# Returns
Parameters related to vapor flux and humidity:
- `DEi`: vapor flux difference [kg/m²s]
- `Eimp_2m`: impervious surface vapor flux at 2m [kg/m²s]
- `Ebare_soil_2m`: bare soil vapor flux at 2m [kg/m²s]
- `Eveg_int_2m`: intercepted water vapor flux at 2m [kg/m²s]
- `Eveg_soil_2m`: soil under vegetation vapor flux at 2m [kg/m²s]
- `TEveg_2m`: vegetation transpiration at 2m [kg/m²s]
- `Ecan_2m`: canyon vapor flux at 2m [kg/m²s]
- `q2m`: specific humidity at 2m [kg/kg]
- `e_T2m`: vapor pressure at 2m [Pa]
- `RH_T2m`: relative humidity at 2m [-]
- `qcan`: canyon specific humidity [kg/kg]
- `e_Tcan`: canyon vapor pressure [Pa]
- `RH_Tcan`: canyon relative humidity [-]
"""
function air_humidity_2m_output(
    q2m::FT,
    T2m::FT,
    Timp::FT,
    Tbare::FT,
    Tveg::FT,
    Tcan::FT,
    qcan::FT,
    rap_can2m::FT,
    rap_can2m_Inv::FT,
    rb_L::FT,
    alp_soil_bare::FT,
    r_soil_bare::FT,
    alp_soil_veg::FT,
    r_soil_veg::FT,
    rs_sun_L::FT,
    rs_shd_L::FT,
    dw_L::FT,
    Fsun_L::FT,
    Fshd_L::FT,
    FractionsGround,
    ParVegGround,
    Eimp::FT,
    Ebare::FT,
    Eveg_int::FT,
    Eveg_pond::FT,
    Eveg_soil::FT,
    TEveg::FT,
    Pre::FT,
    Humidity_ittm,
    fconv::FT,
    MeteoData,
    Gemeotry_m,
    rho_atm::FT,
    Zp1::FT,
    ParCalculation,
) where {FT<:AbstractFloat}
    # Vapor pressure and specific humidity calculations
    esat_Timp = 611 * exp(17.27 * (Timp - 273.16) / (237.3 + (Timp - 273.16)))
    qsat_Timp = (0.622 * esat_Timp) / (Pre - 0.378 * esat_Timp)

    esat_Tbare = 611 * exp(17.27 * (Tbare - 273.16) / (237.3 + (Tbare - 273.16)))
    qsat_Tbare = (0.622 * esat_Tbare) / (Pre - 0.378 * esat_Tbare)

    esat_Tveg = 611 * exp(17.27 * (Tveg - 273.16) / (237.3 + (Tveg - 273.16)))
    qsat_Tveg = (0.622 * esat_Tveg) / (Pre - 0.378 * esat_Tveg)

    esat_Tcan = 611 * exp(17.27 * (Tcan - 273.16) / (237.3 + (Tcan - 273.16)))
    esat_T2m = 611 * exp(17.27 * (T2m - 273.16) / (237.3 + (T2m - 273.16)))

    e_Tcan = qcan * Pre / (0.622 + 0.378 * qcan)
    RH_Tcan = e_Tcan / esat_Tcan

    e_T2m = q2m * Pre / (0.622 + 0.378 * q2m)
    RH_T2m = e_T2m / esat_T2m

    # Enhanced resistance calculation
    if Tcan - MeteoData.Tatm > 0.1
        ra_enhanced = rap_can2m_Inv .* (1 - fconv)
    else
        ra_enhanced = rap_can2m_Inv
    end

    # Turbulent heat fluxes
    Eimp_2m_pot = rho_atm * (qsat_Timp - q2m) ./ rap_can2m
    Eimp_2m = min.(Eimp_2m_pot, Eimp)

    Ebare_soil_2m_pot =
        rho_atm * (alp_soil_bare * qsat_Tbare - q2m) ./ (rap_can2m + r_soil_bare)
    Ebare_soil_2m = min.(Ebare_soil_2m_pot, Ebare)

    Eveg_int_2m_pot =
        rho_atm * (qsat_Tveg - q2m) ./
        (rb_L / ((ParVegGround.LAI + ParVegGround.SAI) * dw_L) + rap_can2m)
    Eveg_int_2m = min.(Eveg_int, Eveg_int_2m_pot)

    Eveg_soil_2m_pot =
        rho_atm * (alp_soil_veg * qsat_Tveg - q2m) ./ (rap_can2m + r_soil_veg)
    Eveg_soil_2m = min.(Eveg_pond + Eveg_soil, Eveg_soil_2m_pot)

    TEveg_sun_2m_pot =
        rho_atm * (qsat_Tveg - q2m) ./ (
            rb_L / ((ParVegGround.LAI) * Fsun_L * (1 - dw_L)) +
            rap_can2m +
            rs_sun_L / ((ParVegGround.LAI) * Fsun_L * (1 - dw_L))
        )
    TEveg_shd_2m_pot =
        rho_atm * (qsat_Tveg - q2m) ./ (
            rb_L / ((ParVegGround.LAI) * Fshd_L * (1 - dw_L)) +
            rap_can2m +
            rs_shd_L / ((ParVegGround.LAI) * Fshd_L * (1 - dw_L))
        )
    TEveg_2m_pot = TEveg_sun_2m_pot + TEveg_shd_2m_pot
    TEveg_2m = min.(TEveg_2m_pot, TEveg)

    Ecan_2m = rho_atm * (q2m - qcan) ./ ra_enhanced

    Vcanyon =
        (
            Gemeotry_m.Width_canyon *
            min(2 * Zp1 / Gemeotry_m.Height_canyon, 1) *
            Gemeotry_m.Height_canyon
        ) / Gemeotry_m.Width_canyon
    dS_E_air = Vcanyon * rho_atm * (q2m - Humidity_ittm.q2m) / ParCalculation.dts

    # Equation set up
    Eimp_2m = FractionsGround.fimp * Eimp_2m
    Ebare_2m = FractionsGround.fbare * Ebare_soil_2m
    Eveg_2m = FractionsGround.fveg * (Eveg_int_2m + Eveg_soil_2m + TEveg_2m)

    DEi = Ecan_2m + dS_E_air - Eimp_2m - Ebare_2m - Eveg_2m

    return DEi,
    Eimp_2m,
    Ebare_soil_2m,
    Eveg_int_2m,
    Eveg_soil_2m,
    TEveg_2m,
    Ecan_2m,
    q2m,
    e_T2m,
    RH_T2m,
    qcan,
    e_Tcan,
    RH_Tcan
end
