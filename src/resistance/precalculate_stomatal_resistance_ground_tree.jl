"""
    precalculate_stomatal_resistance_ground_tree(
        TempVec_ittm::NamedTuple,
        Humidity_ittm::NamedTuple,
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        SoilPotW_ittm::NamedTuple,
        CiCO2Leaf_ittm::NamedTuple,
        MeteoData::NamedTuple,
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
        rb_L::FT,
        rb_H::FT,
        rap_can::FT,
        rap_Htree_In::FT,
    ) where {FT<:AbstractFloat}

Calculate stomatal resistance for ground vegetation and trees.

# Arguments
- `TempVec_ittm`: Temperature vector containing canyon, ground vegetation and tree temperatures
- `Humidity_ittm`: Humidity data structure with canyon specific humidity
- `ParVegGround`: Ground vegetation parameters including LAI, photosynthesis properties
- `SoilPotW_ittm`: Soil water potential data for ground vegetation and trees
- `CiCO2Leaf_ittm`: Leaf CO2 concentration data for previous iteration
- `MeteoData`: Meteorological data including atmospheric pressure and CO2 concentration
- `geometry`: Urban geometry parameters
- `FractionsGround`: Ground surface fractions (vegetated, bare, impervious)
- `PropOpticalGround`: Ground optical properties
- `PropOpticalWall`: Wall optical properties
- `PropOpticalTree`: Tree optical properties
- `ParVegTree`: Tree vegetation parameters
- `SunPosition`: Sun position parameters (solar zenith and azimuth angles)
- `ViewFactor`: View factors between urban surfaces
- `ParWindows`: Window parameters
- `BEM_on`: Building energy model flag
- `rb_L`: Boundary layer resistance for ground vegetation [s/m]
- `rb_H`: Boundary layer resistance for trees [s/m]
- `rap_can`: Within-canyon aerodynamic resistance [s/m]
- `rap_Htree_In`: Aerodynamic resistance within tree canopy [s/m]

# Returns
- `rs_sun_H`: Stomatal resistance for sunlit tree leaves [s/m]
- `rs_shd_H`: Stomatal resistance for shaded tree leaves [s/m]
- `Ci_sun_H`: Internal CO2 concentration for sunlit tree leaves [ppm]
- `Ci_shd_H`: Internal CO2 concentration for shaded tree leaves [ppm]
- `rs_sun_L`: Stomatal resistance for sunlit ground vegetation leaves [s/m]
- `rs_shd_L`: Stomatal resistance for shaded ground vegetation leaves [s/m]
- `Ci_sun_L`: Internal CO2 concentration for sunlit ground vegetation leaves [ppm]
- `Ci_shd_L`: Internal CO2 concentration for shaded ground vegetation leaves [ppm]
"""

function precalculate_stomatal_resistance_ground_tree(
    TempVec_ittm::NamedTuple,
    Humidity_ittm::NamedTuple,
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SoilPotW_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    MeteoData::NamedTuple,
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
    rb_L::FT,
    rb_H::FT,
    rap_can::FT,
    rap_Htree_In::FT,
) where {FT<:AbstractFloat}

    # Calculate shortwave radiation
    _, _, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, _ = Radiation.total_shortwave_absorbed(
        geometry,
        MeteoData.SW_dir,
        MeteoData.SW_diff,
        SunPosition.theta_n,
        SunPosition.theta_Z,
        FractionsGround,
        PropOpticalGround,
        PropOpticalWall,
        PropOpticalTree,
        ParVegTree,
        ViewFactor,
        ParWindows,
        BEM_on,
    )

    # Tree absorbed: conversion from sphere to horizontal projected area
    SWRabs_t_Tree =
        SWRabs_t.Tree * 4 * geometry.radius_tree * π / (4 * geometry.radius_tree)
    SWRabsDir_t_Tree =
        SWRabsDir_t.Tree * 4 * geometry.radius_tree * π / (4 * geometry.radius_tree)
    SWRabsDiff_t_Tree =
        SWRabsDiff_t.Tree * 4 * geometry.radius_tree * π / (4 * geometry.radius_tree)

    SWRdir_abs_tree = SWRabsDir_t_Tree
    SWRdiff_abs_tree = SWRabsDiff_t_Tree
    SWRdir_abs_groundveg = SWRabsDir_t.GroundVeg
    SWRdiff_abs_groundveg = SWRabsDiff_t.GroundVeg

    # Parameter specification
    Tcanyon = TempVec_ittm.TCanyon
    qcanyon = Humidity_ittm.CanyonSpecific
    Tveg = TempVec_ittm.TGroundVeg
    Ttree = TempVec_ittm.TTree
    Pre = MeteoData.Pre
    Catm_O2 = MeteoData.Catm_O2
    Catm_CO2 = MeteoData.Catm_CO2
    fgveg = FractionsGround.fveg
    trees = geometry.trees

    LAI_L = ParVegGround.LAI
    LAI_H = ParVegTree.LAI
    Kopt_H = ParVegTree.Kopt
    Kopt_L = ParVegGround.Kopt
    Knit_H = ParVegTree.Knit
    Knit_L = ParVegGround.Knit
    Psi_sto_50_H = ParVegTree.Psi_sto_50
    Psi_sto_50_L = ParVegGround.Psi_sto_50
    Psi_sto_00_H = ParVegTree.Psi_sto_00
    Psi_sto_00_L = ParVegGround.Psi_sto_00
    CT_H = ParVegTree.CT
    CT_L = ParVegGround.CT
    Vmax_H = ParVegTree.Vmax
    Vmax_L = ParVegGround.Vmax
    DSE_H = ParVegTree.DSE
    DSE_L = ParVegGround.DSE
    Ha_H = ParVegTree.Ha
    Ha_L = ParVegGround.Ha
    FI_H = ParVegTree.FI
    FI_L = ParVegGround.FI
    Do_H = ParVegTree.Do
    Do_L = ParVegGround.Do
    a1_H = ParVegTree.a1
    a1_L = ParVegGround.a1
    go_H = ParVegTree.go
    go_L = ParVegGround.go
    e_rel_H = ParVegTree.e_rel
    e_rel_L = ParVegGround.e_rel
    e_relN_H = ParVegTree.e_relN
    e_relN_L = ParVegGround.e_relN
    gmes_H = ParVegTree.gmes
    gmes_L = ParVegGround.gmes
    rjv_H = ParVegTree.rjv
    rjv_L = ParVegGround.rjv
    mSl_H = ParVegTree.mSl
    mSl_L = ParVegGround.mSl
    Sl_H = ParVegTree.Sl
    Sl_L = ParVegGround.Sl

    Psi_L_tm1 = SoilPotW_ittm.SoilPotWGroundVeg_L
    Psi_H_tm1 = SoilPotW_ittm.SoilPotWGroundTot_H
    Ci_sun_H_tm1 = CiCO2Leaf_ittm.CiCO2LeafTreeSun
    Ci_shd_H_tm1 = CiCO2Leaf_ittm.CiCO2LeafTreeShd
    Ci_sun_L_tm1 = CiCO2Leaf_ittm.CiCO2LeafGroundVegSun
    Ci_shd_L_tm1 = CiCO2Leaf_ittm.CiCO2LeafGroundVegShd

    # Structural parameters
    Cveg = fgveg > 0
    Ctree = trees == 1

    # Vapor pressure saturation and specific humidity at esat
    esat_T_canyon = 611 * exp(17.27 * (Tcanyon - 273.16) / (237.3 + (Tcanyon - 273.16)))
    e_T_canyon = qcanyon * Pre / (0.622 + 0.378 * qcanyon)

    # Parameters for stomata resistance
    Citm1_sun_H = Ctree * Ci_sun_H_tm1
    Citm1_shd_H = Ctree * Ci_shd_H_tm1
    Citm1_sun_L = Cveg * Ci_sun_L_tm1
    Citm1_shd_L = Cveg * Ci_shd_L_tm1

    Ds_canyon = esat_T_canyon - e_T_canyon
    Oa = Catm_O2

    # Partitioning of radiation into sunlit and shaded area for trees
    Fsun_H = Ctree * (1.0 - exp(-Kopt_H * LAI_H)) / (Kopt_H * LAI_H)
    Fsun_H = Fsun_H < 0.01 ? zero(FT) : Fsun_H
    Fsun_H = min(Fsun_H, one(FT))
    Fshd_H = Ctree * (1 - Fsun_H)

    PAR_sun_H = Ctree * (SWRdir_abs_tree + Fsun_H * SWRdiff_abs_tree)
    PAR_shd_H = Ctree * (Fshd_H * SWRdiff_abs_tree)

    # Partitioning for ground vegetation
    Fsun_L = Cveg * (1.0 - exp(-Kopt_L * LAI_L)) / (Kopt_L * LAI_L)
    Fsun_L = Fsun_L < 0.01 ? zero(FT) : Fsun_L
    Fsun_L = min(Fsun_L, one(FT))
    Fshd_L = Cveg * (1 - Fsun_L)

    PAR_sun_L = Cveg * (SWRdir_abs_groundveg + Fsun_L * SWRdiff_abs_groundveg)
    PAR_shd_L = Cveg * (Fshd_L * SWRdiff_abs_groundveg)

    # Calculate stomatal resistances
    Opt_CR = FT(1)  # Numerical tolerance for internal CO2 computation

    # Initialize return values
    rs_sun_H = Inf
    rs_shd_H = Inf
    Ci_sun_H = zero(FT)
    Ci_shd_H = zero(FT)
    rs_sun_L = Inf
    rs_shd_L = Inf
    Ci_sun_L = zero(FT)
    Ci_shd_L = zero(FT)

    if Ctree == 1 && LAI_H > 0
        rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, _, _, _, _ = canopy_resistance_an_evolution(
            PAR_sun_H,
            PAR_shd_H,
            LAI_H,
            Kopt_H,
            Knit_H,
            Fsun_H,
            Fshd_H,
            Citm1_sun_H,
            Citm1_shd_H,
            Catm_CO2,
            rap_Htree_In,
            rb_H,
            Ttree - 273.15,
            Tcanyon - 273.15,
            Pre/100,
            Ds_canyon,
            Psi_H_tm1,
            Psi_sto_50_H,
            Psi_sto_00_H,
            CT_H,
            Vmax_H,
            DSE_H,
            Ha_H,
            FI_H,
            Oa,
            Do_H,
            a1_H,
            go_H,
            e_rel_H,
            e_relN_H,
            gmes_H,
            rjv_H,
            mSl_H,
            Sl_H,
            Opt_CR,
        )
    end

    if Cveg == 1 && LAI_L > 0
        rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L, _, _, _, _ = canopy_resistance_an_evolution(
            PAR_sun_L,
            PAR_shd_L,
            LAI_L,
            Kopt_L,
            Knit_L,
            Fsun_L,
            Fshd_L,
            Citm1_sun_L,
            Citm1_shd_L,
            Catm_CO2,
            rap_can,
            rb_L,
            Tveg - 273.15,
            Tcanyon - 273.15,
            Pre/100,
            Ds_canyon,
            Psi_L_tm1,
            Psi_sto_50_L,
            Psi_sto_00_L,
            CT_L,
            Vmax_L,
            DSE_L,
            Ha_L,
            FI_L,
            Oa,
            Do_L,
            a1_L,
            go_L,
            e_rel_L,
            e_relN_L,
            gmes_L,
            rjv_L,
            mSl_L,
            Sl_L,
            Opt_CR,
        )
    end

    return rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L
end
