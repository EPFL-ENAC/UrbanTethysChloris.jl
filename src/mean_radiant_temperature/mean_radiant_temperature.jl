"""
    mean_radiant_temperature(
        SWRout_t::NamedTuple,
        LWRout_t::NamedTuple,
        MeteoData::NamedTuple,
        ViewFactorPoint::ModelComponents.Parameters.ViewFactorPoint{FT},
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        SunPosition::NamedTuple,
        Person::NamedTuple,
    ) where {FT<:AbstractFloat}

Calculate mean radiant temperature for a person in an urban canyon.

# Arguments
- `SWRout_t`: Outgoing shortwave radiation
- `LWRout_t`: Outgoing longwave radiation
- `MeteoData`: Meteorological data
- `ViewFactorPoint`: View factors from point to surfaces
- `ParTree`: Tree presence parameters
- `ParVegTree`: Tree vegetation parameters
- `Geometry_m`: Urban geometry parameters
- `SunPosition`: Solar position parameters
- `Person`: Person position parameters

# Returns
- `Tmrt::FT`: Mean radiant temperature [°C]
- `BooleanInSun::FT`: Factor indicating if person is in sun [0-1]
- `SWRdir_Person::FT`: Direct shortwave radiation on person [W/m²]
- `SWRdir_in_top::FT`: Direct radiation on top surface [W/m²]
- `SWRdir_in_bottom::FT`: Direct radiation on bottom surface [W/m²]
- `SWRdir_in_east::FT`: Direct radiation on east surface [W/m²]
- `SWRdir_in_south::FT`: Direct radiation on south surface [W/m²]
- `SWRdir_in_west::FT`: Direct radiation on west surface [W/m²]
- `SWRdir_in_north::FT`: Direct radiation on north surface [W/m²]
- `SWRdiff_Person::FT`: Diffuse shortwave radiation on person [W/m²]
- `LWR_Person::FT`: Longwave radiation on person [W/m²]
"""
function mean_radiant_temperature(
    SWRout_t::Radiation.RadiationFluxes{FT},
    LWRout_t::Radiation.RadiationFluxes{FT},
    model::Model{FT},
    ViewFactorPoint::ViewFactorPoint{FT},
) where {FT<:AbstractFloat}
    return mean_radiant_temperature(
        SWRout_t,
        LWRout_t,
        model.forcing.meteorological,
        ViewFactorPoint,
        model.parameters.vegetation.tree,
        model.parameters.urbangeometry,
        model.forcing.sunposition,
        model.parameters.person,
        FT(hour(model.forcing)),
    )
end
function mean_radiant_temperature(
    SWRout_t::Radiation.RadiationFluxes{FT},
    LWRout_t::Radiation.RadiationFluxes{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    ViewFactorPoint::ViewFactorPoint{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    SunPosition::ModelComponents.ForcingInputs.SunPositionInputs{FT,0},
    Person::ModelComponents.Parameters.PersonParameters{FT},
    TimeHr::FT,
) where {FT<:AbstractFloat}

    # Extract parameters
    trees = Geometry_m.trees
    h_can = Geometry_m.hcanyon
    w_can = Geometry_m.wcanyon
    h_tree = Geometry_m.htree
    r_tree = Geometry_m.radius_tree
    d_tree = Geometry_m.distance_tree
    theta_Z = SunPosition.theta_Z
    theta_n = SunPosition.theta_n
    zeta_S = SunPosition.zeta_S
    h_P = Person.PositionPz
    x_P = Person.PositionPx
    SWR_dir = MeteoData.SW_dir
    Wcan = Geometry_m.Width_canyon
    TimeOfMaxSolAlt = SunPosition.TimeOfMaxSolAlt
    # TimeHr = SunPosition.Datam[4]
    SunDSM_MRT = MeteoData.SunDSM_MRT

    # Calculate if person is in shade
    BooleanInSun = person_in_shade(
        trees,
        h_can,
        w_can,
        d_tree,
        h_tree,
        r_tree,
        theta_Z,
        theta_n,
        h_P,
        x_P,
        ParVegTree,
        Wcan,
        TimeOfMaxSolAlt,
        TimeHr,
    )

    # Adjust for DSM shading
    if SunDSM_MRT == 0 && BooleanInSun > 0
        BooleanInSun = zero(FT)
    end

    # Calculate direct and diffuse radiation
    SWRdir_Person, SWRdir_in_top, SWRdir_in_bottom, SWRdir_in_east, SWRdir_in_south, SWRdir_in_west, SWRdir_in_north = swr_dir_person(
        SWR_dir, zeta_S, theta_Z, BooleanInSun
    )

    SWRdiff_Person, LWR_Person = swr_diff_person(
        SWRout_t, LWRout_t, MeteoData, ViewFactorPoint, TimeOfMaxSolAlt, TimeHr
    )

    # Calculate mean radiant temperature
    AbsCoeff = FT(0.7)    # Absorption coefficient for shortwave radiation
    EmCoeff = FT(0.97)    # Emissivity of human body
    bolzm = FT(5.67e-8)   # Stefan-Boltzmann constant [W/m²K⁴]

    # Mean radiant flux density as sum of radiation in three dimensions
    Sstr = AbsCoeff * (SWRdir_Person + SWRdiff_Person) + EmCoeff * LWR_Person

    # Mean radiant temperature from Stefan-Boltzmann law [°C]
    Tmrt = (Sstr / (EmCoeff * bolzm))^(FT(0.25)) - FT(273.25)

    return Tmrt,
    BooleanInSun,
    SWRdir_Person,
    SWRdir_in_top,
    SWRdir_in_bottom,
    SWRdir_in_east,
    SWRdir_in_south,
    SWRdir_in_west,
    SWRdir_in_north,
    SWRdiff_Person,
    LWR_Person
end
function mean_radiant_temperature(
    SWRout_t::NamedTuple,
    LWRout_t::NamedTuple,
    MeteoData::NamedTuple,
    ViewFactorPoint::ViewFactorPoint{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    SunPosition::NamedTuple,
    Person::NamedTuple,
) where {FT<:AbstractFloat}

    # Extract parameters
    trees = Geometry_m.trees
    h_can = Geometry_m.hcanyon
    w_can = Geometry_m.wcanyon
    h_tree = Geometry_m.htree
    r_tree = Geometry_m.radius_tree
    d_tree = Geometry_m.distance_tree
    theta_Z = SunPosition.theta_Z
    theta_n = SunPosition.theta_n
    zeta_S = SunPosition.zeta_S
    h_P = Person.PositionPz
    x_P = Person.PositionPx
    SWR_dir = MeteoData.SW_dir
    Wcan = Geometry_m.Width_canyon
    TimeOfMaxSolAlt = SunPosition.TimeOfMaxSolAlt
    TimeHr = SunPosition.Datam[4]
    SunDSM_MRT = MeteoData.SunDSM_MRT

    # Calculate if person is in shade
    BooleanInSun = person_in_shade(
        trees,
        h_can,
        w_can,
        d_tree,
        h_tree,
        r_tree,
        theta_Z,
        theta_n,
        h_P,
        x_P,
        ParVegTree,
        Wcan,
        TimeOfMaxSolAlt,
        TimeHr,
    )

    # Adjust for DSM shading
    if SunDSM_MRT == 0 && BooleanInSun > 0
        BooleanInSun = zero(FT)
    end

    # Calculate direct and diffuse radiation
    SWRdir_Person, SWRdir_in_top, SWRdir_in_bottom, SWRdir_in_east, SWRdir_in_south, SWRdir_in_west, SWRdir_in_north = swr_dir_person(
        SWR_dir, zeta_S, theta_Z, BooleanInSun
    )

    SWRdiff_Person, LWR_Person = swr_diff_person(
        SWRout_t, LWRout_t, MeteoData, ViewFactorPoint, TimeOfMaxSolAlt, TimeHr
    )

    # Calculate mean radiant temperature
    AbsCoeff = FT(0.7)    # Absorption coefficient for shortwave radiation
    EmCoeff = FT(0.97)    # Emissivity of human body
    bolzm = FT(5.67e-8)   # Stefan-Boltzmann constant [W/m²K⁴]

    # Mean radiant flux density as sum of radiation in three dimensions
    Sstr = AbsCoeff * (SWRdir_Person + SWRdiff_Person) + EmCoeff * LWR_Person

    # Mean radiant temperature from Stefan-Boltzmann law [°C]
    Tmrt = (Sstr / (EmCoeff * bolzm))^(FT(0.25)) - FT(273.25)

    return Tmrt,
    BooleanInSun,
    SWRdir_Person,
    SWRdir_in_top,
    SWRdir_in_bottom,
    SWRdir_in_east,
    SWRdir_in_south,
    SWRdir_in_west,
    SWRdir_in_north,
    SWRdiff_Person,
    LWR_Person
end
