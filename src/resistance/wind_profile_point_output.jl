"""
    wind_profile_point_output(
        Zp::FT,
        Gemeotry_m::NamedTuple,
        ParVegTree::NamedTuple,
        ParTree::NamedTuple,
        MeteoData::NamedTuple,
        FractionsGround::NamedTuple,
        ParVegGround::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate wind speed at a specific height in urban canyon.

# Arguments
- `Zp`: Height of interest [m]
- `Gemeotry_m`: Canyon geometry parameters
- `ParVegTree`: Tree vegetation parameters
- `ParTree`: Tree presence parameters
- `MeteoData`: Meteorological data
- `FractionsGround`: Ground fraction parameters
- `ParVegGround`: Ground vegetation parameters

# Returns
- `u_Zp`: Wind speed at height Zp [m/s]
"""
function wind_profile_point_output!(model::Model{FT}) where {FT<:AbstractFloat}
    u_ZPerson = wind_profile_point_output(
        model.parameters.person.HeightWind,
        model.parameters.urbangeometry,
        model.parameters.vegetation.tree,
        model.forcing.meteorological,
        model.parameters.surfacefractions.ground,
        model.parameters.vegetation.ground,
    )

    model.variables.environmentalconditions.wind.u_ZPerson = u_ZPerson

    return nothing
end

function wind_profile_point_output(
    Zp::FT,
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
) where {FT<:AbstractFloat}
    # Extract parameters
    Hcan = Gemeotry_m.Height_canyon
    Wcan = Gemeotry_m.Width_canyon
    Wroof = Gemeotry_m.Width_roof
    Htree = Gemeotry_m.Height_tree
    R_tree = Gemeotry_m.Radius_tree
    Hcan_max = Gemeotry_m.Hcan_max
    Hcan_std = Gemeotry_m.Hcan_std

    Kopt = ParVegTree.Kopt
    LAI_t = ParVegTree.LAI
    trees = Gemeotry_m.trees

    Zatm = MeteoData.Zatm
    uatm = MeteoData.Uatm

    fgveg = FractionsGround.fveg
    fgbare = FractionsGround.fbare
    fgimp = FractionsGround.fimp

    hc_L = ParVegGround.hc
    Zref_und = FT(1.5)

    # Calculate presence indicators
    Cimp = fgimp > 0
    Cbare = fgbare > 0
    Cveg = fgveg > 0
    Ctree = trees

    # Get ground roughness parameters
    _, _, zom_ground, _, _, _, _, _, _, _, _, _ = urban_roughness(
        Ctree * Htree, Cveg * hc_L, Cbare, Cimp, false
    )

    # Calculate wind profile
    _, _, _, u_Zp, _, _, _ = wind_profile_canyon(
        Hcan,
        Htree,
        R_tree,
        Wcan,
        Wroof,
        Kopt,
        LAI_t,
        Zatm,
        uatm,
        Zp,
        trees,
        Zref_und,
        zom_ground,
        Hcan_max,
        Hcan_std,
    )

    return u_Zp
end
function wind_profile_point_output(
    Zp::FT,
    Gemeotry_m::NamedTuple,
    ParVegTree::NamedTuple,
    ParTree::NamedTuple,
    MeteoData::NamedTuple,
    FractionsGround::NamedTuple,
    ParVegGround::NamedTuple,
) where {FT<:AbstractFloat}
    # Extract parameters
    Hcan = Gemeotry_m.Height_canyon
    Wcan = Gemeotry_m.Width_canyon
    Wroof = Gemeotry_m.Width_roof
    Htree = Gemeotry_m.Height_tree
    R_tree = Gemeotry_m.Radius_tree
    Hcan_max = Gemeotry_m.Hcan_max
    Hcan_std = Gemeotry_m.Hcan_std

    Kopt = ParVegTree.Kopt
    LAI_t = ParVegTree.LAI
    trees = ParTree.trees

    Zatm = MeteoData.Zatm
    uatm = MeteoData.Uatm

    fgveg = FractionsGround.fveg
    fgbare = FractionsGround.fbare
    fgimp = FractionsGround.fimp

    hc_L = ParVegGround.hc
    Zref_und = FT(1.5)

    # Calculate presence indicators
    Cimp = fgimp > 0
    Cbare = fgbare > 0
    Cveg = fgveg > 0
    Ctree = trees

    # Get ground roughness parameters
    _, _, zom_ground, _, _, _, _, _, _, _, _, _ = urban_roughness(
        Ctree * Htree, Cveg * hc_L, Cbare, Cimp, false
    )

    # Calculate wind profile
    _, _, _, u_Zp, _, _, _ = wind_profile_canyon(
        Hcan,
        Htree,
        R_tree,
        Wcan,
        Wroof,
        Kopt,
        LAI_t,
        Zatm,
        uatm,
        Zp,
        trees,
        Zref_und,
        zom_ground,
        Hcan_max,
        Hcan_std,
    )

    return u_Zp
end
