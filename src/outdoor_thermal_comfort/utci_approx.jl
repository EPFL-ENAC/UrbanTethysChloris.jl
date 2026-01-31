"""
    utci_approx(
        Ta::FT,
        RH::FT,
        Tmrt::FT,
        va::FT
    ) where {FT<:AbstractFloat}

Calculate Universal Thermal Climate Index (UTCI) using polynomial approximation.

# Arguments
- `Ta`: Air temperature [°C]
- `RH`: Relative humidity [%]
- `Tmrt`: Mean radiant temperature [°C]
- `va`: Wind speed at 10m height [m/s]

# Returns
- `UTCI::FT`: Universal Thermal Climate Index [°C]
"""
function utci_approx!(model::Model{FT}) where {FT<:AbstractFloat}
    T2m = model.variables.humidity.Results2m.T2m
    RH_T2m = model.variables.humidity.Results2m.RH_T2m
    Tmrt = model.variables.temperature.mrt.Tmrt
    u_ZPerson = model.variables.environmentalconditions.wind.u_ZPerson

    UTCI = utci_approx(T2m - FT(273.15), RH_T2m * 100, Tmrt, u_ZPerson)

    model.variables.temperature.thermalcomfort.UTCI = UTCI

    return nothing
end

function utci_approx(Ta::FT, RH::FT, Tmrt::FT, va::FT) where {FT<:AbstractFloat}
    esat = saturation_vapor_pressure(Ta)
    ehPa = esat * RH / 100.0

    D_Tmrt = Tmrt - Ta
    Pa = ehPa / 10.0

    va_square = va * va
    va_cube = va_square * va
    va_quad = va_square * va_square

    Ta_square = Ta * Ta
    Ta_cube = Ta_square * Ta
    Ta_quad = Ta_square * Ta_square
    Ta_pent = Ta_quad * Ta
    Ta_hex = Ta_pent * Ta

    D_Tmrt_square = D_Tmrt * D_Tmrt
    D_Tmrt_cube = D_Tmrt_square * D_Tmrt
    D_Tmrt_quad = D_Tmrt_square * D_Tmrt_square

    Pa_square = Pa * Pa
    Pa_cube = Pa_square * Pa

    UTCI_approx =
        Ta +
        FT(6.07562052e-01) +
        FT(-2.27712343e-02) * Ta +
        FT(8.06470249e-04) * Ta_square +
        FT(-1.54271372e-04) * Ta_cube +
        FT(-3.24651735e-06) * Ta_quad +
        FT(7.32602852e-08) * Ta_pent +
        FT(1.35959073e-09) * Ta_hex +
        FT(-2.25836520e+00) * va +
        FT(8.80326035e-02) * Ta * va +
        FT(2.16844454e-03) * Ta_square * va +
        FT(-1.53347087e-05) * Ta_cube * va +
        FT(-5.72983704e-07) * Ta_quad * va +
        FT(-2.55090145e-09) * Ta_pent * va +
        FT(-7.51269505e-01) * va_square +
        FT(-4.08350271e-03) * Ta * va_square +
        FT(-5.21670675e-05) * Ta_square * va_square +
        FT(1.94544667e-06) * Ta_cube * va_square +
        FT(1.14099531e-08) * Ta_quad * va_square +
        FT(1.58137256e-01) * va_cube +
        FT(-6.57263143e-05) * Ta * va_cube +
        FT(2.22697524e-07) * Ta_square * va_cube +
        FT(-4.16117031e-08) * Ta_cube * va_cube +
        FT(-1.27762753e-02) * va_quad +
        FT(9.66891875e-06) * Ta * va_quad +
        FT(2.52785852e-09) * Ta_square * va_quad +
        FT(4.56306672e-04) * va_quad * va +
        FT(-1.74202546e-07) * Ta * va_quad * va +
        FT(-5.91491269e-06) * va_quad * va_square +
        FT(3.98374029e-01) * D_Tmrt +
        FT(1.83945314e-04) * Ta * D_Tmrt +
        FT(-1.73754510e-04) * Ta_square * D_Tmrt +
        FT(-7.60781159e-07) * Ta_cube * D_Tmrt +
        FT(3.77830287e-08) * Ta_quad * D_Tmrt +
        FT(5.43079673e-10) * Ta_pent * D_Tmrt +
        FT(-2.00518269e-02) * va * D_Tmrt +
        FT(8.92859837e-04) * Ta * va * D_Tmrt +
        FT(3.45433048e-06) * Ta_square * va * D_Tmrt +
        FT(-3.77925774e-07) * Ta_cube * va * D_Tmrt +
        FT(-1.69699377e-09) * Ta_quad * va * D_Tmrt +
        FT(1.69992415e-04) * va_square * D_Tmrt +
        FT(-4.99204314e-05) * Ta * va_square * D_Tmrt +
        FT(2.47417178e-07) * Ta_square * va_square * D_Tmrt +
        FT(1.07596466e-08) * Ta_cube * va_square * D_Tmrt +
        FT(8.49242932e-05) * va_cube * D_Tmrt +
        FT(1.35191328e-06) * Ta * va_cube * D_Tmrt +
        FT(-6.21531254e-09) * Ta_square * va_cube * D_Tmrt +
        FT(-4.99410301e-06) * va_quad * D_Tmrt +
        FT(-1.89489258e-08) * Ta * va_quad * D_Tmrt +
        FT(8.15300114e-08) * va_quad * va * D_Tmrt +
        FT(7.55043090e-04) * D_Tmrt_square +
        FT(-5.65095215e-05) * Ta * D_Tmrt_square +
        FT(-4.52166564e-07) * Ta_square * D_Tmrt_square +
        FT(2.46688878e-08) * Ta_cube * D_Tmrt_square +
        FT(2.42674348e-10) * Ta_quad * D_Tmrt_square +
        FT(1.54547250e-04) * va * D_Tmrt_square +
        FT(5.24110970e-06) * Ta * va * D_Tmrt_square +
        FT(-8.75874982e-08) * Ta_square * va * D_Tmrt_square +
        FT(-1.50743064e-09) * Ta_cube * va * D_Tmrt_square +
        FT(-1.56236307e-05) * va_square * D_Tmrt_square +
        FT(-1.33895614e-07) * Ta * va_square * D_Tmrt_square +
        FT(2.49709824e-09) * Ta_square * va_square * D_Tmrt_square +
        FT(6.51711721e-07) * va_cube * D_Tmrt_square +
        FT(1.94960053e-09) * Ta * va_cube * D_Tmrt_square +
        FT(-1.00361113e-08) * va_quad * D_Tmrt_square +
        FT(-1.21206673e-05) * D_Tmrt_cube +
        FT(-2.18203660e-07) * Ta * D_Tmrt_cube +
        FT(7.51269482e-09) * Ta_square * D_Tmrt_cube +
        FT(9.79063848e-11) * Ta_cube * D_Tmrt_cube +
        FT(1.25006734e-06) * va * D_Tmrt_cube +
        FT(-1.81584736e-09) * Ta * va * D_Tmrt_cube +
        FT(-3.52197671e-10) * Ta_square * va * D_Tmrt_cube +
        FT(-3.36514630e-08) * va_square * D_Tmrt_cube +
        FT(1.35908359e-10) * Ta * va_square * D_Tmrt_cube +
        FT(4.17032620e-10) * va_cube * D_Tmrt_cube +
        FT(-1.30369025e-09) * D_Tmrt_quad +
        FT(4.13908461e-10) * Ta * D_Tmrt_quad +
        FT(9.22652254e-12) * Ta_square * D_Tmrt_quad +
        FT(-5.08220384e-09) * va * D_Tmrt_quad +
        FT(-2.24730961e-11) * Ta * va * D_Tmrt_quad +
        FT(1.17139133e-10) * va_square * D_Tmrt_quad +
        FT(6.62154879e-10) * D_Tmrt_quad * D_Tmrt +
        FT(4.03863260e-13) * Ta * D_Tmrt_quad * D_Tmrt +
        FT(1.95087203e-12) * va * D_Tmrt_quad * D_Tmrt +
        FT(-4.73602469e-12) * D_Tmrt_cube * D_Tmrt_quad +
        FT(5.12733497e+00) * Pa +
        FT(-3.12788561e-01) * Ta * Pa +
        FT(-1.96701861e-02) * Ta_square * Pa +
        FT(9.99690870e-04) * Ta_cube * Pa +
        FT(9.51738512e-06) * Ta_quad * Pa +
        FT(-4.66426341e-07) * Ta_pent * Pa +
        FT(5.48050612e-01) * va * Pa +
        FT(-3.30552823e-03) * Ta * va * Pa +
        FT(-1.64119440e-03) * Ta_square * va * Pa +
        FT(-5.16670694e-06) * Ta_cube * va * Pa +
        FT(9.52692432e-07) * Ta_quad * va * Pa +
        FT(-4.29223622e-02) * va_square * Pa +
        FT(5.00845667e-03) * Ta * va_square * Pa +
        FT(1.00601257e-06) * Ta_square * va_square * Pa +
        FT(-1.81748644e-06) * Ta_cube * va_square * Pa +
        FT(-1.25813502e-03) * va_cube * Pa +
        FT(-1.79330391e-04) * Ta * va_cube * Pa +
        FT(2.34994441e-06) * Ta_square * va_cube * Pa +
        FT(1.29735808e-04) * va_quad * Pa +
        FT(1.29064870e-06) * Ta * va_quad * Pa +
        FT(-2.28558686e-06) * va_quad * va * Pa +
        FT(-3.69476348e-02) * D_Tmrt * Pa +
        FT(1.62325322e-03) * Ta * D_Tmrt * Pa +
        FT(-3.14279680e-05) * Ta_square * D_Tmrt * Pa +
        FT(2.59835559e-06) * Ta_cube * D_Tmrt * Pa +
        FT(-4.77136523e-08) * Ta_quad * D_Tmrt * Pa +
        FT(8.64203390e-03) * va * D_Tmrt * Pa +
        FT(-6.87405181e-04) * Ta * va * D_Tmrt * Pa +
        FT(-9.13863872e-06) * Ta_square * va * D_Tmrt * Pa +
        FT(5.15916806e-07) * Ta_cube * va * D_Tmrt * Pa +
        FT(-3.59217476e-05) * va_square * D_Tmrt * Pa +
        FT(3.28696511e-05) * Ta * va_square * D_Tmrt * Pa +
        FT(-7.10542454e-07) * Ta_square * va_square * D_Tmrt * Pa +
        FT(-1.24382300e-05) * va_cube * D_Tmrt * Pa +
        FT(-7.38584400e-09) * Ta * va_cube * D_Tmrt * Pa +
        FT(2.20609296e-07) * va_quad * D_Tmrt * Pa +
        FT(-7.32469180e-04) * D_Tmrt_square * Pa +
        FT(-1.87381964e-05) * Ta * D_Tmrt_square * Pa +
        FT(4.80925239e-06) * Ta_square * D_Tmrt_square * Pa +
        FT(-8.75492040e-08) * Ta_cube * D_Tmrt_square * Pa +
        FT(2.77862930e-05) * va * D_Tmrt_square * Pa +
        FT(-5.06004592e-06) * Ta * va * D_Tmrt_square * Pa +
        FT(1.14325367e-07) * Ta_square * va * D_Tmrt_square * Pa +
        FT(2.53016723e-06) * va_square * D_Tmrt_square * Pa +
        FT(-1.72857035e-08) * Ta * va_square * D_Tmrt_square * Pa +
        FT(-3.95079398e-08) * va_cube * D_Tmrt_square * Pa +
        FT(-3.59413173e-07) * D_Tmrt_cube * Pa +
        FT(7.04388046e-07) * Ta * D_Tmrt_cube * Pa +
        FT(-1.89309167e-08) * Ta_square * D_Tmrt_cube * Pa +
        FT(-4.79768731e-07) * va * D_Tmrt_cube * Pa +
        FT(7.96079978e-09) * Ta * va * D_Tmrt_cube * Pa +
        FT(1.62897058e-09) * va_square * D_Tmrt_cube * Pa +
        FT(3.94367674e-08) * D_Tmrt_quad * Pa +
        FT(-1.18566247e-09) * Ta * D_Tmrt_quad * Pa +
        FT(3.34678041e-10) * va * D_Tmrt_quad * Pa +
        FT(-1.15606447e-10) * D_Tmrt_quad * D_Tmrt * Pa +
        FT(-2.80626406e+00) * Pa_square +
        FT(5.48712484e-01) * Ta * Pa_square +
        FT(-3.99428410e-03) * Ta_square * Pa_square +
        FT(-9.54009191e-04) * Ta_cube * Pa_square +
        FT(1.93090978e-05) * Ta_quad * Pa_square +
        FT(-3.08806365e-01) * va * Pa_square +
        FT(1.16952364e-02) * Ta * va * Pa_square +
        FT(4.95271903e-04) * Ta_square * va * Pa_square +
        FT(-1.90710882e-05) * Ta_cube * va * Pa_square +
        FT(2.10787756e-03) * va_square * Pa_square +
        FT(-6.98445738e-04) * Ta * va_square * Pa_square +
        FT(2.30109073e-05) * Ta_square * va_square * Pa_square +
        FT(4.17856590e-04) * va_cube * Pa_square +
        FT(-1.27043871e-05) * Ta * va_cube * Pa_square +
        FT(-3.04620472e-06) * va_quad * Pa_square +
        FT(5.14507424e-02) * D_Tmrt * Pa_square +
        FT(-4.32510997e-03) * Ta * D_Tmrt * Pa_square +
        FT(8.99281156e-05) * Ta_square * D_Tmrt * Pa_square +
        FT(-7.14663943e-07) * Ta_cube * D_Tmrt * Pa_square +
        FT(-2.66016305e-04) * va * D_Tmrt * Pa_square +
        FT(2.63789586e-04) * Ta * va * D_Tmrt * Pa_square +
        FT(-7.01199003e-06) * Ta_square * va * D_Tmrt * Pa_square +
        FT(-1.06823306e-04) * va_square * D_Tmrt * Pa_square +
        FT(3.61341136e-06) * Ta * va_square * D_Tmrt * Pa_square +
        FT(2.29748967e-07) * va_cube * D_Tmrt * Pa_square +
        FT(3.04788893e-04) * D_Tmrt_square * Pa_square +
        FT(-6.42070836e-05) * Ta * D_Tmrt_square * Pa_square +
        FT(1.16257971e-06) * Ta_square * D_Tmrt_square * Pa_square +
        FT(7.68023384e-06) * va * D_Tmrt_square * Pa_square +
        FT(-5.47446896e-07) * Ta * va * D_Tmrt_square * Pa_square +
        FT(-3.59937910e-08) * va_square * D_Tmrt_square * Pa_square +
        FT(-4.36497725e-06) * D_Tmrt_cube * Pa_square +
        FT(1.68737969e-07) * Ta * D_Tmrt_cube * Pa_square +
        FT(2.67489271e-08) * va * D_Tmrt_cube * Pa_square +
        FT(3.23926897e-09) * D_Tmrt_quad * Pa_square +
        FT(-3.53874123e-02) * Pa_cube +
        FT(-2.21201190e-01) * Ta * Pa_cube +
        FT(1.55126038e-02) * Ta_square * Pa_cube +
        FT(-2.63917279e-04) * Ta_cube * Pa_cube +
        FT(4.53433455e-02) * va * Pa_cube +
        FT(-4.32943862e-03) * Ta * va * Pa_cube +
        FT(1.45389826e-04) * Ta_square * va * Pa_cube +
        FT(2.17508610e-04) * va_square * Pa_cube +
        FT(-6.66724702e-05) * Ta * va_square * Pa_cube +
        FT(3.33217140e-05) * va_cube * Pa_cube +
        FT(-2.26921615e-03) * D_Tmrt * Pa_cube +
        FT(3.80261982e-04) * Ta * D_Tmrt * Pa_cube +
        FT(-5.45314314e-09) * Ta_square * D_Tmrt * Pa_cube +
        FT(-7.96355448e-04) * va * D_Tmrt * Pa_cube +
        FT(2.53458034e-05) * Ta * va * D_Tmrt * Pa_cube +
        FT(-6.31223658e-06) * va_square * D_Tmrt * Pa_cube +
        FT(3.02122035e-04) * D_Tmrt_square * Pa_cube +
        FT(-4.77403547e-06) * Ta * D_Tmrt_square * Pa_cube +
        FT(1.73825715e-06) * va * D_Tmrt_square * Pa_cube +
        FT(-4.09087898e-07) * D_Tmrt_cube * Pa_cube +
        FT(6.14155345e-01) * Pa_square * Pa_square +
        FT(-6.16755931e-02) * Ta * Pa_square * Pa_square +
        FT(1.33374846e-03) * Ta_square * Pa_square * Pa_square +
        FT(3.55375387e-03) * va * Pa_square * Pa_square +
        FT(-5.13027851e-04) * Ta * va * Pa_square * Pa_square +
        FT(1.02449757e-04) * va_square * Pa_square * Pa_square +
        FT(-1.48526421e-03) * D_Tmrt * Pa_square * Pa_square +
        FT(-4.11469183e-05) * Ta * D_Tmrt * Pa_square * Pa_square +
        FT(-6.80434415e-06) * va * D_Tmrt * Pa_square * Pa_square +
        FT(-9.77675906e-06) * D_Tmrt_square * Pa_square * Pa_square +
        FT(8.82773108e-02) * Pa_cube * Pa_square +
        FT(-3.01859306e-03) * Ta * Pa_cube * Pa_square +
        FT(1.04452989e-03) * va * Pa_cube * Pa_square +
        FT(2.47090539e-04) * D_Tmrt * Pa_cube * Pa_square +
        FT(1.48348065e-03) * Pa_cube * Pa_cube

    return UTCI_approx
end

"""
    saturation_vapor_pressure(ta::FT) where {FT<:AbstractFloat}

Calculate saturation vapor pressure over water in hPa for input air temperature in Celsius.

# Arguments
- `ta`: Air temperature [°C]

# Returns
- `esat::FT`: Saturation vapor pressure [hPa]
"""
function saturation_vapor_pressure(ta::FT) where {FT<:AbstractFloat}
    # Constants from Hardy (1998)
    gV = [
        FT(-2.8365744e3),
        FT(-6.028076559e3),
        FT(1.954263612e1),
        FT(-2.737830188e-2),
        FT(1.6261698e-5),
        FT(7.0229056e-10),
        FT(-1.8680009e-13),
        FT(2.7150305),
    ]

    tk = ta + FT(273.15)  # Convert to Kelvin
    es = gV[8] * log(tk)

    # Calculate polynomial sum
    for i in 1:7
        es += gV[i] * tk^(i - 3)
    end

    # Convert Pa to hPa
    esat = exp(es) * FT(0.01)

    return esat
end
