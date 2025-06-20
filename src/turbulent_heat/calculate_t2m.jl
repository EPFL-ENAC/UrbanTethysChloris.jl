"""
    calculate_t2m(
        Timp::FT,
        Tbare::FT,
        Tveg::FT,
        Twsun::FT,
        Twshade::FT,
        Tcan::FT,
        Zp1::FT,
        rap_can2m::FT,
        rap_can2m_Inv::FT,
        rb_L::FT,
        RES_w1::FT,
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        Gemeotry_m::NamedTuple,
        geometry::NamedTuple,
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        TempVec_ittm,
        cp_atm::FT,
        rho_atm::FT,
        ParCalculation::NamedTuple,
        fconv::FT,
        MeteoData::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate the air temperature at 2m level and compute sensible heat fluxes.

# Arguments
- `Timp`: impervious ground temperature [K]
- `Tbare`: bare ground temperature [K]
- `Tveg`: vegetated ground temperature [K]
- `Twsun`: sunlit wall temperature [K]
- `Twshade`: shaded wall temperature [K]
- `Tcan`: canyon air temperature [K]
- `Zp1`: reference height [m]
- `rap_can2m`: resistance between canyon and 2m height [s/m]
- `rap_can2m_Inv`: inverse resistance between canyon and 2m height [m/s]
- `rb_L`: leaf boundary layer resistance [s/m]
- `RES_w1`: wall resistance [s/m]
- `FractionsGround`: ground cover fractions
- `Gemeotry_m`: urban geometry parameters
- `geometry`: urban geometry parameters
- `ParVegGround`: vegetation parameters
- `TempVec_ittm`: temperature variables from previous timestep
- `cp_atm`: specific heat capacity of air [J/kg/K]
- `rho_atm`: air density [kg/m³]
- `ParCalculation`: calculation parameters
- `fconv`: convection factor [-]
- `MeteoData`: meteorological data

# Returns
- `T2m::FT`: air temperature at 2m [K]
- `DHi::FT`: sensible heat flux difference [W/m²]
- `Himp_2m::FT`: impervious surface sensible heat flux at 2m [W/m²]
- `Hbare_2m::FT`: bare soil sensible heat flux at 2m [W/m²]
- `Hveg_2m::FT`: vegetated surface sensible heat flux at 2m [W/m²]
- `Hwsun_2m::FT`: sunlit wall sensible heat flux at 2m [W/m²]
- `Hwshade_2m::FT`: shaded wall sensible heat flux at 2m [W/m²]
- `Hcan_2m::FT`: canyon sensible heat flux at 2m [W/m²]
"""
function calculate_t2m(
    Timp::FT,
    Tbare::FT,
    Tveg::FT,
    Twsun::FT,
    Twshade::FT,
    Tcan::FT,
    Zp1::FT,
    rap_can2m::FT,
    rap_can2m_Inv::FT,
    rb_L::FT,
    RES_w1::FT,
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    Gemeotry_m::NamedTuple,
    geometry::NamedTuple,
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    TempVec_ittm,
    cp_atm::FT,
    rho_atm::FT,
    ParCalculation::NamedTuple,
    fconv::FT,
    MeteoData::NamedTuple,
) where {FT<:AbstractFloat}
    Vcanyon =
        (
            Gemeotry_m.Width_canyon *
            min(2 * Zp1 / Gemeotry_m.Height_canyon, 1) *
            Gemeotry_m.Height_canyon
        ) / Gemeotry_m.Width_canyon

    if Tcan - MeteoData.Tatm > 0.1
        ra_enhanced = rap_can2m_Inv * (1 - fconv)
    else
        ra_enhanced = rap_can2m_Inv
    end

    # Resistance calculations
    Rimp_H = FractionsGround.fimp / rap_can2m
    Rbare_H = FractionsGround.fbare / rap_can2m
    Rveg_H =
        FractionsGround.fveg /
        (rb_L / (2 * (ParVegGround.LAI + ParVegGround.SAI)) + rap_can2m)
    Rwsun_H = min(2 * Zp1 / Gemeotry_m.Height_canyon, 1) * geometry.hcanyon / RES_w1
    Rwshd_H = min(2 * Zp1 / Gemeotry_m.Height_canyon, 1) * geometry.hcanyon / RES_w1
    Rcan2m_H = 1 / ra_enhanced
    RdS_H = Vcanyon / ParCalculation.dts

    # Temperature at 2m calculation
    T2m =
        (
            Rcan2m_H * Tcan +
            RdS_H * TempVec_ittm.T2m +
            Rimp_H * Timp +
            Rbare_H * Tbare +
            Rveg_H * Tveg +
            Rwsun_H * Twsun +
            Rwshd_H * Twshade
        ) / (Rcan2m_H + RdS_H + Rimp_H + Rbare_H + Rveg_H + Rwsun_H + Rwshd_H)

    # Heat flux calculations
    Himp_2m = FractionsGround.fimp * cp_atm * rho_atm * ((Timp - T2m) / rap_can2m)
    Hbare_2m = FractionsGround.fbare * cp_atm * rho_atm * ((Tbare - T2m) / rap_can2m)
    Hveg_2m =
        FractionsGround.fveg *
        cp_atm *
        rho_atm *
        ((Tveg - T2m) / (rb_L / (2 * (ParVegGround.LAI + ParVegGround.SAI)) + rap_can2m))
    Hwsun_2m =
        min(2 * Zp1 / Gemeotry_m.Height_canyon, 1) *
        geometry.hcanyon *
        cp_atm *
        rho_atm *
        ((Twsun - T2m) / RES_w1)
    Hwshade_2m =
        min(2 * Zp1 / Gemeotry_m.Height_canyon, 1) *
        geometry.hcanyon *
        cp_atm *
        rho_atm *
        ((Twshade - T2m) / RES_w1)
    Hcan_2m = cp_atm * rho_atm * (T2m - Tcan) / rap_can2m_Inv
    dS_H_air = Vcanyon * cp_atm * rho_atm * (T2m - TempVec_ittm.T2m) / ParCalculation.dts

    DHi = Hcan_2m + dS_H_air - Himp_2m - Hbare_2m - Hveg_2m - Hwsun_2m - Hwshade_2m

    return T2m, DHi, Himp_2m, Hbare_2m, Hveg_2m, Hwsun_2m, Hwshade_2m, Hcan_2m
end
