"""
    soil_resistance(
        Ts::FT,
        Pre::FT,
        Ws::FT,
        ea::FT,
        q_runon::FT,
        O::FT,
        Ks::FT,
        Osat::FT,
        Ohy::FT,
        L::FT,
        Pe::FT,
        O33::FT,
        alpVG::FT,
        nVG::FT,
        SPAR::Int
    ) where {FT<:AbstractFloat}

Calculate soil resistance parameters.

# Arguments
- `Ts`: Soil temperature [°C]
- `Pre`: Atmospheric pressure [hPa]
- `Ws`: Wind speed [m/s]
- `ea`: Vapor pressure [Pa]
- `q_runon`: Runon [mm/dth]
- `O`: Water content [-]
- `Ks`: Saturated hydraulic conductivity [mm/h]
- `Osat`: Saturated water content [-]
- `Ohy`: Residual water content [-]
- `L`: Pore size distribution index [-]
- `Pe`: Air entry pressure [kPa]
- `O33`: Water content at -33 kPa [-]
- `alpVG`: van Genuchten alpha parameter [1/mm]
- `nVG`: van Genuchten n parameter [-]
- `SPAR`: Soil parameter type (1-VanGenuchten, 2-Saxton-Rawls)

# Returns
- `r_soil`: Soil resistance [s/m]
- `b_soil`: Beta factor [-]
- `alp_soil`: Relative humidity [-]
"""

function soil_resistance(
    Ts, Pre, Ws, ea, q_runon, O, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, SPAR
)
    Ts_k = Ts + 273.15  # Soil Temperature [K]

    # Constants
    row = 1000.0    # [kg/m^3] water density
    g = 9.81        # [m/s^2] gravity acceleration
    Rd = 461.5      # [J/kg K] water vapor gas constant
    Da = (2.11e-5) * ((Ts_k / 273.15)^1.94) * (Pre * 100 / 101325)  # [m^2/s] vapor molecular diffusivity
    esat = 611.0 * exp.(17.27 .* Ts ./ (237.3 .+ Ts))  # [Pa]

    Ko, Po = Soil.conductivity_suction(SPAR, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, O)
    Po = max.(Po, 0)  # [mm]
    alp_soil = exp.(-Po * g / (1000 * Rd * Ts_k))  # [-]

    # From Haghighi et al 2013
    Psz = (11.12 * nVG^3.286) * 1e-6  # [m] 40-200 um  Size of the pores
    dm = 2.26e-3 / sqrt(Ws)  # Shahraeeni et al 2012 [m] Boundary Layer Thickness

    ANSW = 4
    if ANSW == 1
        # r_soil = exp(8.206 - 4.255*(O-Ohy)/(Osat-Ohy))
        # r_soil = O <= Ohy ? Inf : r_soil
        # b_soil = 1

    elseif ANSW == 2
        # b_soil = 0.25 * (1 - cos(O/O33*π))^2
        # b_soil = O >= O33 ? 1 : b_soil
        # b_soil = O <= Ohy ? 0 : b_soil
        # r_soil = 0

    elseif ANSW == 3
        # # Omitting case 3 as it uses undefined variables (ro, qa, ra, etc.)
        # r_soil = 0
        # b_soil = 1

    elseif ANSW == 4
        gammap = (alp_soil * esat - ea) / (row * Rd * Ts_k)  # [-]
        if gammap < 0
            r_soil = 0.0
        else
            rsv = gammap / (4 * Ko / (1000 * 3600))  # Internal soil viscous resistance [s/m]
            f_O = (2 / π) * (sqrt(π / (4 * O)) - 1) / sqrt(4 * O)  # [-]
            rvbl = (dm + Psz * f_O) / Da  # [s/m] viscous boundary layer resistance
            r_soil = rvbl + rsv
        end
        r_soil = O <= Ohy ? Inf : r_soil
        b_soil = 1.0

    elseif ANSW == 5
        # # From Lehmann et al 2018
        # mVG = 1 - 1/nVG
        # hc = -1/alpVG * ((nVG-1)/nVG)^((1-2*nVG)/nVG)  # [mm]
        # # SPAR != 1 && @warn "ERROR: OPTION FOR SOIL RESISTANCE REQUIRES van-Genuchten soil hydraulic parameters"

        # gammap = (alp_soil*esat-ea)/(row*Rd*Ts_k)  # [-]
        # E0 = Da/dm*gammap*1000/3600  # [mm/h]
        # Se = 1/((1+abs(alpVG*hc)^nVG)^mVG)
        # Se = min(Se, 1)
        # Khc = Ks*((Se)^(0.5))*(1-(1-(Se)^(1/mVG))^mVG)^2  # [mm/h]
        # dhdz = 1 + E0/(4*Khc)  # [-]

        # if gammap < 0
        #     r_soil = 0
        # else
        #     rsv = gammap/(4*Ko*dhdz/(1000*3600))
        #     f_O = (2/π)*(sqrt(π/(4*O))-1)/sqrt(4*O)
        #     rvbl = (dm + Psz*f_O)/Da
        #     r_soil = rvbl + rsv
        # end
        # r_soil = O <= Ohy ? Inf : r_soil
        # b_soil = 1
    end

    if q_runon > 0
        r_soil = 0.0
        alp_soil = 1.0
        b_soil = 1.0
    end

    return r_soil, b_soil, alp_soil
end
