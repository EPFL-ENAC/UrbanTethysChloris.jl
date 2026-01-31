"""
    post_calculate_soil_moisture_change(
        OwaterInitial::ModelComponents.ModelVariables.Owater{FT},
        Owater::ModelComponents.ModelVariables.Owater{FT},
        ParSoilRoof::ModelComponents.Parameters.SoilParameters{FT},
        ParSoilGround::ModelComponents.Parameters.SoilParameters{FT},
        FractionsRoof::ModelComponents.Parameters.SurfaceFractionsParameters{FT},
        FractionsGround::ModelComponents.Parameters.SurfaceFractionsParameters{FT},
        geometry::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    ) where {FT<:AbstractFloat}

Post-calculates soil water volume changes for roof, canyon, and urban areas.

# Arguments
- `OwaterInitial`: Initial soil water content at the start of simulation.
- `Owater`: Current soil water content.
- `ParSoilRoof`: Soil parameters for roof.
- `ParSoilGround`: Soil parameters for ground/canyon.
- `FractionsRoof`: Surface fractions for roof.
- `FractionsGround`: Surface fractions for ground/canyon.
- `geometry`: Urban geometry parameters.

# Returns
- `dVdtRoofCalc`: Calculated change in soil water volume for roof (mm/time step).
- `dVdtCanCalc`: Calculated change in soil water volume for canyon (mm/time step).
- `dVdtUrbCalc`: Calculated change in soil water volume for urban area (mm/time step).
"""
function post_calculate_soil_moisture_change(
    OwaterInitial::ModelComponents.ModelVariables.Owater{FT},
    results::Dict{String,Any},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    geometry::ModelComponents.Parameters.UrbanGeometryParameters{FT},
) where {FT<:AbstractFloat}
    function _soil_parameters(
        soil::ModelComponents.Parameters.VegetatedSoilParameters{FT}
    ) where {FT<:AbstractFloat}
        _, dz, _, _, Ohy = Soil.soil_parameters_total(
            soil.Pcla,
            soil.Psan,
            soil.Porg,
            soil.Kfc,
            soil.Phy,
            soil.SPAR,
            soil.Kbot,
            1,
            1,
            zeros(FT, 1),
            zeros(FT, 1),
            zeros(FT, 1),
            zeros(FT, 1),
            zeros(FT, 1),
            zeros(FT, 1),
            soil.Zs,
        )

        return dz, Ohy
    end

    ParSoilRoofdz, ParSoilRoofOhy = _soil_parameters(ParSoilRoof)
    ParSoilGrounddz, ParSoilGroundOhy = _soil_parameters(ParSoilGround)

    # @infiltrate
    # Initial soil water setting in beginning of simulation (time step 1)
    Vinit_Rveg = (OwaterInitial.OwRoofSoilVeg .- ParSoilRoofOhy) .* ParSoilRoofdz
    Vinit_Gimp = (OwaterInitial.OwGroundSoilImp .- ParSoilGroundOhy) .* ParSoilGrounddz
    Vinit_Gbare = (OwaterInitial.OwGroundSoilBare .- ParSoilGroundOhy) .* ParSoilGrounddz
    Vinit_Gveg = (OwaterInitial.OwGroundSoilVeg .- ParSoilGroundOhy) .* ParSoilGrounddz

    # Water volume in each soil column
    VRveg = (results["OwRoofSoilVeg"] .- ParSoilRoofOhy) .* ParSoilRoofdz
    VGimp = (results["OwGroundSoilImp"] .- ParSoilGroundOhy) .* ParSoilGrounddz
    VGbare = (results["OwGroundSoilBare"] .- ParSoilGroundOhy) .* ParSoilGrounddz
    VGveg = (results["OwGroundSoilVeg"] .- ParSoilGroundOhy) .* ParSoilGrounddz

    # Total water volume in soil (canyon & roof)
    VinitRoof = FractionsRoof.fveg * NaNMath.sum(Vinit_Rveg)
    VinitCan =
        FractionsGround.fimp * NaNMath.sum(Vinit_Gimp) +
        FractionsGround.fbare * NaNMath.sum(Vinit_Gbare) +
        FractionsGround.fveg * NaNMath.sum(Vinit_Gveg)

    VRoof = FractionsRoof.fveg * sum(VRveg)
    VCan =
        FractionsGround.fimp * sum(VGimp) +
        FractionsGround.fbare * sum(VGbare) +
        FractionsGround.fveg * sum(VGveg)

    # Change in soil volume in soil column
    dVdtRoofCalc = VRoof - VinitRoof
    dVdtCanCalc = VCan - VinitCan
    dVdtUrbCalc = geometry.wcanyon_norm * dVdtCanCalc + geometry.wroof_norm * dVdtRoofCalc

    return dVdtRoofCalc, dVdtCanCalc, dVdtUrbCalc
end
