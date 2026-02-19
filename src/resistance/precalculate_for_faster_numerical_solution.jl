"""
    precalculate_for_faster_numerical_solution(
        model::Model{FT},
        TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
        Humidity_ittm::ModelComponents.ModelVariables.Humidity{FT},
        SoilPotW_ittm::ModelComponents.ModelVariables.SoilPotW{FT},
        CiCO2Leaf_ittm::ModelComponents.ModelVariables.CiCO2Leaf{FT},
        RES_ittm::ModelComponents.ModelVariables.Resistance{FT},
        ittn::Int,
        ittm::Int,
        ViewFactor::RayTracing.ViewFactor{FT},
        BEM_on::Bool,
    ) where {FT<:AbstractFloat}

Calculate enhancement factor and precalculate stomatal resistances for faster numerical solution.

# Arguments
- `model::Model{FT}`: The model instance containing variables and parameters
- `TempVec_ittm`: Temperature vector containing surface temperatures
- `Humidity_ittm`: Humidity data structure
- `SoilPotW_ittm`: Soil water potential data
- `CiCO2Leaf_ittm`: Leaf CO2 concentration data
- `RES_ittm`: Resistance values from previous iterations
- `ittn::Int`: Current time step index
- `ittm::Int`: Current sensitivity index
- `ViewFactor`: View factors between urban surfaces
- `BEM_on`: Building energy model flag

# Returns
- `fconv`: Enhancement factor for convective resistance [-]
- `rsRoofPreCalc`: StomatalResistancePreCalc with roof vegetation parameters
- `rsGroundPreCalc`: StomatalResistancePreCalc with ground vegetation parameters
- `rsTreePreCalc`: StomatalResistancePreCalc with tree vegetation parameters
"""
function precalculate_for_faster_numerical_solution(
    model::Model{FT},
    model_ittm::ModelIttm{FT,MR,MG},
    ittn::Int,
    ittm::Int,
    ViewFactor::RayTracing.ViewFactor{FT},
    options::AbstractModelOptions,
) where {FT<:AbstractFloat,MR,MG}
    return precalculate_for_faster_numerical_solution(
        ittn,
        ittm,
        model_ittm.tempvec,
        model_ittm.humidity,
        model.parameters.vegetation.ground,
        model_ittm.soilpotw,
        model_ittm.cico2leaf,
        model.forcing.meteorological,
        model.parameters.urbangeometry,
        model.parameters.surfacefractions.ground,
        model.parameters.optical.ground,
        model.parameters.optical.wall,
        model.parameters.optical.tree,
        model.parameters.vegetation.tree,
        model.forcing.sunposition,
        ViewFactor,
        model.parameters.building_energy.windows,
        options.BEM_on,
        model.parameters.vegetation.roof,
        model.parameters.optical.roof,
        model.parameters.surfacefractions.roof,
        model_ittm.resistance,
        options.OPT_Obhukov,
    )
end

# TODO: remove references to ittm
function precalculate_for_faster_numerical_solution(
    ittn::Int,
    ittm::Int,
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    Humidity_ittm::ModelComponents.ModelVariables.Humidity{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SoilPotW_ittm::ModelComponents.ModelVariables.SoilPotW{FT},
    CiCO2Leaf_ittm::ModelComponents.ModelVariables.CiCO2Leaf{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    geometry::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    PropOpticalTree::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SunPosition::ModelComponents.ForcingInputs.SunPositionInputs{Float64,0},
    ViewFactor::RayTracing.ViewFactor{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    RES_ittm::ModelComponents.ModelVariables.Resistance{FT},
    strategy::AbstractZeroFindingStrategies=SimpleBrentStrategy(FT),
) where {FT<:AbstractFloat}

    # Calculate enhancement factor based on ra_original of previous time step
    if ittn == 1
        fconv = zero(FT)
    else
        if TempVec_ittm.TCanyon - TempVec_ittm.Tatm > 0.1
            hPBL = FT(1000)  # Planetary boundary layer height [m]

            dcan, zomcan, _, _, _, _, _ = wind_profile_canyon(
                geometry.Height_canyon,
                geometry.Height_tree,
                geometry.Radius_tree,
                geometry.Width_canyon,
                geometry.Width_roof,
                ParVegTree.Kopt,
                ParVegTree.LAI,
                MeteoData.Zatm,
                MeteoData.Uatm,
                geometry.Height_canyon,
                geometry.trees,
                FT(1.5),
                FT(0.01),
                geometry.Hcan_max,
                geometry.Hcan_std,
            )

            zom_town = zomcan                  # Momentum roughness length of canyon, calculated according to McDonald 1998
            zoh_town = zom_town / FT(10)       # Heat roughness length of canyon

            fconv, _, _, _ = enhancement_factor_ra_pleim(
                RES_ittm.raCanyontoAtmOrig,
                zom_town,
                zoh_town,
                dcan,
                MeteoData.Zatm,
                MeteoData.Uatm,
                hPBL,
                strategy,
            )
        else
            fconv = zero(FT)
        end
    end

    # Precalculate stomatal resistance for roof
    if ParVegRoof.LAI > 0 && FractionsRoof.fveg > 0
        if ittn == 1
            ra = FT(100)
            rb = FT(50)
            # Stomatal resistance, roughly 200-300 during the day and ca 3000 during the night
            rs_sun, rs_shd, Ci_sun, Ci_shd = precalculate_stomatal_resistance_roof(
                TempVec_ittm,
                MeteoData,
                Humidity_ittm,
                ParVegRoof,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                PropOpticalRoof,
                ra,
                rb,
            )
        else
            rs_sun, rs_shd, Ci_sun, Ci_shd = precalculate_stomatal_resistance_roof(
                TempVec_ittm,
                MeteoData,
                Humidity_ittm,
                ParVegRoof,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                PropOpticalRoof,
                RES_ittm.raRooftoAtm,
                RES_ittm.rb_LRoof,
            )
        end
    else
        rs_sun = Inf
        rs_shd = Inf
        Ci_sun = zero(FT)
        Ci_shd = zero(FT)
    end

    rsRoofPreCalc = StomatalResistancePreCalc{FT}(;
        rs_sun=rs_sun, rs_shd=rs_shd, Ci_sun=Ci_sun, Ci_shd=Ci_shd
    )

    # Precalculate stomatal resistance for tree and for ground vegetation in canyon
    if ittn == 1
        rb_LGround = FT(50)
        rb_HGround = FT(50)
        rap_can = FT(100)
        rap_Htree_In = FT(100)

        rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L = precalculate_stomatal_resistance_ground_tree(
            TempVec_ittm,
            Humidity_ittm,
            ParVegGround,
            SoilPotW_ittm,
            CiCO2Leaf_ittm,
            MeteoData,
            geometry,
            FractionsGround,
            PropOpticalGround,
            PropOpticalWall,
            PropOpticalTree,
            ParVegTree,
            SunPosition,
            ViewFactor,
            ParWindows,
            BEM_on,
            rb_LGround,
            rb_HGround,
            rap_can,
            rap_Htree_In,
        )
    else
        rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L = precalculate_stomatal_resistance_ground_tree(
            TempVec_ittm,
            Humidity_ittm,
            ParVegGround,
            SoilPotW_ittm,
            CiCO2Leaf_ittm,
            MeteoData,
            geometry,
            FractionsGround,
            PropOpticalGround,
            PropOpticalWall,
            PropOpticalTree,
            ParVegTree,
            SunPosition,
            ViewFactor,
            ParWindows,
            BEM_on,
            RES_ittm.rb_LGround,
            RES_ittm.rb_HGround,
            RES_ittm.rap_can,
            RES_ittm.rap_Htree_In,
        )
    end

    rsGroundPreCalc = StomatalResistancePreCalc{FT}(;
        rs_sun=rs_sun_L, rs_shd=rs_shd_L, Ci_sun=Ci_sun_L, Ci_shd=Ci_shd_L
    )

    rsTreePreCalc = StomatalResistancePreCalc{FT}(;
        rs_sun=rs_sun_H, rs_shd=rs_shd_H, Ci_sun=Ci_sun_H, Ci_shd=Ci_shd_H
    )

    return fconv, rsRoofPreCalc, rsGroundPreCalc, rsTreePreCalc
end

function precalculate_for_faster_numerical_solution(
    ittn::Int,
    ittm::Int,
    TempVec_ittm::NamedTuple,
    Humidity_ittm::NamedTuple,
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SoilPotW_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    MeteoData::NamedTuple,
    HumidityAtm::NamedTuple,
    geometry::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    PropOpticalGround::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    PropOpticalWall::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    PropOpticalTree::ModelComponents.Parameters.SimpleOpticalProperties{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SunPosition::NamedTuple,
    ViewFactor::RayTracing.ViewFactor{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    RES::NamedTuple,
    strategy::AbstractZeroFindingStrategies=SimpleBrentStrategy(FT),
) where {FT<:AbstractFloat}

    # Calculate enhancement factor based on ra_original of previous time step
    if ittn == 1
        fconv = zero(FT)
    else
        if TempVec_ittm.TCanyon - TempVec_ittm.Tatm > 0.1
            hPBL = FT(1000)  # Planetary boundary layer height [m]

            dcan, zomcan, _, _, _, _, _ = wind_profile_canyon(
                geometry.Height_canyon,
                geometry.Height_tree,
                geometry.Radius_tree,
                geometry.Width_canyon,
                geometry.Width_roof,
                ParVegTree.Kopt,
                ParVegTree.LAI,
                MeteoData.Zatm,
                MeteoData.Uatm,
                geometry.Height_canyon,
                geometry.trees,
                FT(1.5),
                FT(0.01),
                geometry.Hcan_max,
                geometry.Hcan_std,
            )

            zom_town = zomcan                  # Momentum roughness length of canyon, calculated according to McDonald 1998
            zoh_town = zom_town / FT(10)       # Heat roughness length of canyon

            fconv, _, _, _ = enhancement_factor_ra_pleim(
                RES.raCanyontoAtmOrig[ittn - 1, 1, ittm],
                zom_town,
                zoh_town,
                dcan,
                MeteoData.Zatm,
                MeteoData.Uatm,
                hPBL,
                strategy,
            )
        else
            fconv = zero(FT)
        end
    end

    # Precalculate stomatal resistance for roof
    if ParVegRoof.LAI > 0 && FractionsRoof.fveg > 0
        if ittn == 1
            ra = FT(100)
            rb = FT(50)
            # Stomatal resistance, roughly 200-300 during the day and ca 3000 during the night
            rs_sun, rs_shd, Ci_sun, Ci_shd = precalculate_stomatal_resistance_roof(
                TempVec_ittm,
                MeteoData,
                HumidityAtm,
                ParVegRoof,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                PropOpticalRoof,
                ra,
                rb,
            )
        else
            rs_sun, rs_shd, Ci_sun, Ci_shd = precalculate_stomatal_resistance_roof(
                TempVec_ittm,
                MeteoData,
                HumidityAtm,
                ParVegRoof,
                SoilPotW_ittm,
                CiCO2Leaf_ittm,
                PropOpticalRoof,
                RES.raRooftoAtm[ittn - 1, 1, ittm],
                RES.rb_LRoof[ittn - 1, 1, ittm],
            )
        end
    else
        rs_sun = Inf
        rs_shd = Inf
        Ci_sun = zero(FT)
        Ci_shd = zero(FT)
    end

    rsRoofPreCalc = StomatalResistancePreCalc{FT}(;
        rs_sun=rs_sun, rs_shd=rs_shd, Ci_sun=Ci_sun, Ci_shd=Ci_shd
    )

    # Precalculate stomatal resistance for tree and for ground vegetation in canyon
    if ittn == 1
        rb_LGround = FT(50)
        rb_HGround = FT(50)
        rap_can = FT(100)
        rap_Htree_In = FT(100)

        rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L = precalculate_stomatal_resistance_ground_tree(
            TempVec_ittm,
            Humidity_ittm,
            ParVegGround,
            SoilPotW_ittm,
            CiCO2Leaf_ittm,
            MeteoData,
            geometry,
            FractionsGround,
            PropOpticalGround,
            PropOpticalWall,
            PropOpticalTree,
            ParVegTree,
            SunPosition,
            ViewFactor,
            ParWindows,
            BEM_on,
            rb_LGround,
            rb_HGround,
            rap_can,
            rap_Htree_In,
        )
    else
        rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L = precalculate_stomatal_resistance_ground_tree(
            TempVec_ittm,
            Humidity_ittm,
            ParVegGround,
            SoilPotW_ittm,
            CiCO2Leaf_ittm,
            MeteoData,
            geometry,
            FractionsGround,
            PropOpticalGround,
            PropOpticalWall,
            PropOpticalTree,
            ParVegTree,
            SunPosition,
            ViewFactor,
            ParWindows,
            BEM_on,
            RES.rb_LGround[ittn - 1, 1, ittm],
            RES.rb_HGround[ittn - 1, 1, ittm],
            RES.rap_can[ittn - 1, 1, ittm],
            RES.rap_Htree_In[ittn - 1, 1, ittm],
        )
    end

    rsGroundPreCalc = StomatalResistancePreCalc{FT}(;
        rs_sun=rs_sun_L, rs_shd=rs_shd_L, Ci_sun=Ci_sun_L, Ci_shd=Ci_shd_L
    )

    rsTreePreCalc = StomatalResistancePreCalc{FT}(;
        rs_sun=rs_sun_H, rs_shd=rs_shd_H, Ci_sun=Ci_sun_H, Ci_shd=Ci_shd_H
    )

    return fconv, rsRoofPreCalc, rsGroundPreCalc, rsTreePreCalc
end
