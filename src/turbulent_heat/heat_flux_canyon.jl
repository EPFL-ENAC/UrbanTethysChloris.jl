"""
    heat_flux_canyon(
        TemperatureC::AbstractVector{FT},
        Gemeotry_m,
        MeteoData,
        ParVegTree,
        ParTree,
        fconvPreCalc::FT,
        fconv::FT
    ) where {FT<:AbstractFloat}

Calculate sensible and latent heat fluxes in the urban canyon.

# Arguments
- `TemperatureC`: Temperature vector [K]
- `Gemeotry_m`: Urban geometry parameters
- `MeteoData`: Meteorological data
- `ParVegTree`: Tree vegetation parameters
- `fconvPreCalc`: Pre-calculated convection flag [-]
- `fconv`: Convection factor [-]
"""
function heat_flux_canyon(
    TemperatureC::AbstractVector{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    MeteoData::NamedTuple,
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    fconvPreCalc::FT,
    fconv::FT,
) where {FT<:AbstractFloat}
    T_canyon = TemperatureC[9]
    q_canyon = TemperatureC[10]

    # Extract geometry parameters
    H = Gemeotry_m.Height_canyon
    W = Gemeotry_m.Width_canyon
    Wroof = Gemeotry_m.Width_roof
    Htree = Gemeotry_m.Height_tree
    R_tree = Gemeotry_m.Radius_tree
    Hcan_max = Gemeotry_m.Hcan_max
    Hcan_std = Gemeotry_m.Hcan_std

    # Extract vegetation parameters
    Kopt = ParVegTree.Kopt
    LAI_t = ParVegTree.LAI
    trees = Gemeotry_m.trees

    # Extract meteorological data
    Zatm = MeteoData.Zatm
    Tatm = MeteoData.Tatm
    Uatm = MeteoData.Uatm
    Pre = MeteoData.Pre
    q_atm = MeteoData.q_atm
    ea = MeteoData.ea

    # Calculate thermodynamic properties
    cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    rho_atm = (Pre / (287.04 * Tatm)) * (1 - (ea / Pre) * (1 - 0.622))
    L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    # Vapor pressure and humidity calculations
    esat_T_canyon = 611 * exp(17.27 * (T_canyon - 273.16) / (237.3 + (T_canyon - 273.16)))
    qsat_T_canyon = (0.622 * esat_T_canyon) / (Pre - 0.378 * esat_T_canyon)
    e_T_canyon = q_canyon * Pre / (0.622 + 0.378 * q_canyon)
    rel_hum_canyon = e_T_canyon / esat_T_canyon

    # Calculate structural parameters and wind profile
    dcan, zomcan, _, _, _, _, _ = wind_profile_canyon(
        H,
        Htree,
        R_tree,
        W,
        Wroof,
        Kopt,
        LAI_t,
        Zatm,
        Uatm,
        H,
        trees,
        FT(1.5),
        FT(0.01),
        Hcan_max,
        Hcan_std,
    )

    zom_town = zomcan
    zoh_town = zom_town / 10

    # Calculate aerodynamic resistance
    ra = aerodynamic_resistance(
        Tatm - 273.15,
        T_canyon - 273.15,
        Pre / 100,
        Zatm,
        dcan,
        zom_town,
        zoh_town,
        Uatm,
        ea,
        e_T_canyon,
    )

    ra_orig = ra

    # Include enhancement term for aerodynamic resistance
    if fconvPreCalc == 1
        if T_canyon - Tatm > 0.1
            ra_enhanced = ra * (1 - fconv)
        else
            ra_enhanced = ra
        end
    else
        if T_canyon - Tatm > 0.1
            hPBL = FT(1000)
            fconv, ra_enhanced, _, _ = enhancement_factor_ra_pleim(
                ra, zom_town, zoh_town, dcan, Zatm, Uatm, hPBL
            )
            ra_enhanced = ra * (1 - fconv)
        else
            ra_enhanced = ra
        end
    end

    ra = ra_enhanced

    # Calculate heat fluxes
    ra_canyon = ra
    Hcanyon = cp_atm * rho_atm * (T_canyon - Tatm) / ra
    Ecanyon = rho_atm * (q_canyon - q_atm) / ra
    LEcanyon = L_heat * Ecanyon

    # Create humidity output structure
    HumidityCan = (
        CanyonRelative=rel_hum_canyon,
        CanyonSpecific=q_canyon,
        CanyonVapourPre=e_T_canyon,
        CanyonRelativeSat=one(FT),
        CanyonSpecificSat=qsat_T_canyon,
        CanyonVapourPreSat=esat_T_canyon,
    )

    return Hcanyon, LEcanyon, ra_canyon, ra_orig, fconv, HumidityCan
end
