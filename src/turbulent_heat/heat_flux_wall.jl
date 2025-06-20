"""
    heat_flux_wall(
        TemperatureC::AbstractVector{FT},
        Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        MeteoData::NamedTuple,
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT}
    ) where {FT<:AbstractFloat}

Calculate sensible and latent heat fluxes for sunlit and shaded walls.

# Arguments
- `TemperatureC`: Temperature vector [K]
- `Gemeotry_m`: Urban geometry parameters
- `MeteoData`: Meteorological data
- `ParVegTree`: Tree vegetation parameters
- `ParVegGround`: Ground vegetation parameters
- `FractionsGround`: Ground cover fractions

# Returns
- `Hwsun::FT`: Sensible heat flux from sunlit wall [W/m²]
- `Hwshade::FT`: Sensible heat flux from shaded wall [W/m²]
- `Ewsun::FT`: Latent heat flux from sunlit wall [kg/m²s]
- `Ewshade::FT`: Latent heat flux from shaded wall [kg/m²s]
- `LEwsun::FT`: Latent heat from sunlit wall [W/m²]
- `LEwshade::FT`: Latent heat from shaded wall [W/m²]
- `RES_w1::FT`: Wall resistance at level 1 [s/m]
- `RES_w2::FT`: Wall resistance at level 2 [s/m]
- `rap_Zp1_In::FT`: In-canyon aerodynamic resistance at level 1 [s/m]
- `rap_Zp2_In::FT`: In-canyon aerodynamic resistance at level 2 [s/m]
- `Hwsun1::FT`: Sensible heat flux from sunlit wall at level 1 [W/m²]
- `Hwshade1::FT`: Sensible heat flux from shaded wall at level 1 [W/m²]
- `Hwsun2::FT`: Sensible heat flux from sunlit wall at level 2 [W/m²]
- `Hwshade2::FT`: Sensible heat flux from shaded wall at level 2 [W/m²]
- `cp_atm::FT`: Specific heat capacity of air [J/kg/K]
- `rho_atm::FT`: Air density [kg/m³]
- `L_heat::FT`: Latent heat of vaporization [J/kg]
- `Zp1::FT`: Height of level 1 [m]
- `Zp2::FT`: Height of level 2 [m]
- `rap_Zp1::FT`: Aerodynamic resistance at level 1 [s/m]
- `rap_Zp2::FT`: Aerodynamic resistance at level 2 [s/m]
"""
function heat_flux_wall(
    TemperatureC::AbstractVector{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    MeteoData::NamedTuple,
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
) where {FT<:AbstractFloat}
    # Extract temperatures
    Timp = TemperatureC[1]
    Tbare = TemperatureC[2]
    Tveg = TemperatureC[3]
    Ttree = TemperatureC[6]
    Twsun = TemperatureC[4]
    Twshade = TemperatureC[5]
    T_canyon = TemperatureC[9]
    qcanyon = TemperatureC[10]

    # Extract geometry parameters
    H = Gemeotry_m.Height_canyon
    W = Gemeotry_m.Width_canyon
    Wroof = Gemeotry_m.Width_roof
    Htree = Gemeotry_m.Height_tree
    R_tree = Gemeotry_m.Radius_tree
    Hcan_max = Gemeotry_m.Hcan_max
    Hcan_std = Gemeotry_m.Hcan_std
    rad_tree = R_tree / W

    # Extract vegetation parameters
    Kopt = ParVegTree.Kopt
    LAI_t = ParVegTree.LAI
    trees = Gemeotry_m.trees
    hc_L = ParVegGround.hc

    # Extract meteorological data
    Zatm = MeteoData.Zatm
    Tatm = MeteoData.Tatm
    Uatm = MeteoData.Uatm
    Pre = MeteoData.Pre
    ea = MeteoData.ea

    # Calculate thermodynamic properties
    cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    rho_atm = (Pre / (287.04 * Tatm)) * (1 - (ea / Pre) * (1 - 0.622))
    L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    # Vapor pressure calculation
    e_T_canyon = qcanyon * Pre / (0.622 + 0.378 * qcanyon)

    # Average surface temperature calculation
    Ctree = FT(trees)
    Tsurf =
        (
            FractionsGround.fveg * Tveg +
            FractionsGround.fbare * Tbare +
            FractionsGround.fimp * Timp +
            Ctree * Ttree * (4 * rad_tree) +
            H/W * Twsun +
            H/W * Twshade
        ) / (
            FractionsGround.fveg +
            FractionsGround.fbare +
            FractionsGround.fimp +
            (Ctree * 4 * rad_tree) +
            2 * H/W
        )

    # Calculate roughness parameters
    roughness = urban_roughness(
        zero(FT),
        FT(FractionsGround.fveg > 0) * hc_L,
        FractionsGround.fbare > 0,
        FractionsGround.fimp > 0,
        false,
    )
    zom_ground = roughness[3]

    # Calculate wind profile parameters
    wind_params = wind_profile_canyon(
        H,
        Htree,
        R_tree,
        W,
        Wroof,
        Kopt,
        LAI_t,
        Zatm,
        Uatm,
        FT(2),
        trees,
        FT(1.5),
        zom_ground,
        Hcan_max,
        Hcan_std,
    )
    dcan = wind_params[1]
    zomcan = wind_params[2]
    RoughnessParameter = wind_params[7]

    # Define calculation heights
    Zp1 = FT(2)
    Zp2 = 2 * Zp1 + (H - 2 * Zp1) / 2
    wcan = zero(FT)

    # Calculate in-canyon aerodynamic resistances
    canyon_res = in_canyon_aerodynamic_resistance(
        Uatm,
        Zatm,
        T_canyon-273.15,
        Tsurf-273.15,
        Hcan_max,
        H,
        dcan,
        zomcan,
        FT(1.5),
        zom_ground,
        Zp1,
        Zp2,
        FT(2),
        Pre,
        e_T_canyon,
        RoughnessParameter,
    )
    rap_Zp1 = canyon_res[2]
    rap_Zp1_In = canyon_res[3]
    rap_Zp2 = canyon_res[4]
    rap_Zp2_In = canyon_res[5]
    u_Zp1 = canyon_res[9]
    u_Zp2 = canyon_res[10]

    # Calculate wall resistances
    RES_w1 = cp_atm * rho_atm / (11.8 + 4.2 * sqrt(u_Zp1^2 + wcan^2))
    RES_w2 = cp_atm * rho_atm / (11.8 + 4.2 * sqrt(u_Zp2^2 + wcan^2))

    # Calculate heat fluxes for each level
    Hwsun1 = cp_atm * rho_atm * (Twsun - T_canyon) / (RES_w1 + rap_Zp1_In)
    Hwshade1 = cp_atm * rho_atm * (Twshade - T_canyon) / (RES_w1 + rap_Zp1_In)

    Hwsun2 = cp_atm * rho_atm * (Twsun - T_canyon) / (RES_w2 + rap_Zp2_In)
    Hwshade2 = cp_atm * rho_atm * (Twshade - T_canyon) / (RES_w2 + rap_Zp2_In)

    # Combine fluxes
    Hwsun = min(2 * Zp1/H, 1) * Hwsun1 + max((H - 2 * Zp1)/H, 0) * Hwsun2
    Hwshade = min(2 * Zp1/H, 1) * Hwshade1 + max((H - 2 * Zp1)/H, 0) * Hwshade2

    # No latent heat for walls
    Ewsun = zero(FT)
    Ewshade = zero(FT)
    LEwsun = L_heat * Ewsun
    LEwshade = L_heat * Ewshade

    return Hwsun,
    Hwshade,
    Ewsun,
    Ewshade,
    LEwsun,
    LEwshade,
    RES_w1,
    RES_w2,
    rap_Zp1_In,
    rap_Zp2_In,
    Hwsun1,
    Hwshade1,
    Hwsun2,
    Hwshade2,
    cp_atm,
    rho_atm,
    L_heat,
    Zp1,
    Zp2,
    rap_Zp1,
    rap_Zp2
end
