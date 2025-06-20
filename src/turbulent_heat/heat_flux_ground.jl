"""
    heat_flux_ground(
        TemperatureC::AbstractVector{FT},
        TempVec_ittm::NamedTuple,
        MeteoData::NamedTuple,
        Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        SoilPotW_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        Vwater_ittm::NamedTuple,
        ExWater_ittm::NamedTuple,
        Int_ittm::NamedTuple,
        CiCO2Leaf_ittm::NamedTuple,
        ParInterceptionTree::NamedTuple,
        ParCalculation::NamedTuple,
        SWRdir_abs_tree::FT,
        SWRdiff_abs_tree::FT,
        SWRdir_abs_groundveg::FT,
        SWRdiff_abs_groundveg::FT,
        RESPreCalc::Bool,
        rsGroundPreCalc::NamedTuple,
        rsTreePreCalc::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate sensible and latent heat fluxes for ground surfaces.

# Arguments
- `TemperatureC`: Temperature vector [K]
- `TempVec_ittm`: Temperature from previous timestep
- `MeteoData`: Meteorological data
- `Gemeotry_m`: Urban geometry parameters
- `FractionsGround`: Ground surface fractions
- `ParVegGround`: Ground vegetation parameters
- `ParVegTree`: Tree vegetation parameters
- `ParSoilGround`: Ground soil parameters
- `SoilPotW_ittm`: Soil water potential from previous timestep
- `Owater_ittm`: Soil water content from previous timestep
- `Vwater_ittm`: Soil water volume from previous timestep
- `ExWater_ittm`: Extractable water from previous timestep
- `Int_ittm`: Interception from previous timestep
- `CiCO2Leaf_ittm`: Leaf CO2 concentration from previous timestep
- `ParInterceptionTree`: Tree interception parameters
- `ParCalculation`: Calculation parameters
- `SWRdir_abs_tree`: Direct shortwave radiation absorbed by trees [W/m²]
- `SWRdiff_abs_tree`: Diffuse shortwave radiation absorbed by trees [W/m²]
- `SWRdir_abs_groundveg`: Direct shortwave radiation absorbed by ground vegetation [W/m²]
- `SWRdiff_abs_groundveg`: Diffuse shortwave radiation absorbed by ground vegetation [W/m²]
- `RESPreCalc`: Use pre-calculated resistances
- `rsGroundPreCalc`: Pre-calculated ground resistance parameters
- `rsTreePreCalc`: Pre-calculated tree resistance parameters

# Returns
A NamedTuple containing:
- `Himp::FT`: Sensible heat flux from impervious ground [W/m²]
- `Hbare::FT`: Sensible heat flux from bare ground [W/m²]
- `Hveg::FT`: Sensible heat flux from vegetated ground [W/m²]
- `Htree::FT`: Sensible heat flux from trees [W/m²]
- `Eimp::FT`: Water vapor flux from impervious ground [kg/m²s]
- `Ebare_pond::FT`: Ponded water vapor flux from bare ground [kg/m²s]
- `Ebare_soil::FT`: Soil water vapor flux from bare ground [kg/m²s]
- `Eveg_int::FT`: Intercepted water vapor flux from vegetation [kg/m²s]
- `Eveg_pond::FT`: Ponded water vapor flux from vegetated ground [kg/m²s]
- `Eveg_soil::FT`: Soil water vapor flux from vegetated ground [kg/m²s]
- `TEveg::FT`: Transpiration flux from ground vegetation [kg/m²s]
- `Etree_int::FT`: Intercepted water vapor flux from trees [kg/m²s]
- `TEtree::FT`: Transpiration flux from trees [kg/m²s]
- And many more components including latent heat fluxes and resistances...
"""
function heat_flux_ground(
    TemperatureC::AbstractVector{FT},
    TempVec_ittm::NamedTuple,
    MeteoData::NamedTuple,
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    SoilPotW_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    Vwater_ittm::NamedTuple,
    ExWater_ittm::NamedTuple,
    Int_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    ParInterceptionTree::NamedTuple,
    ParCalculation::NamedTuple,
    SWRdir_abs_tree::FT,
    SWRdiff_abs_tree::FT,
    SWRdir_abs_groundveg::FT,
    SWRdiff_abs_groundveg::FT,
    RESPreCalc::Bool,
    rsGroundPreCalc::NamedTuple,
    rsTreePreCalc::NamedTuple,
) where {FT<:AbstractFloat}
    # Extract temperatures
    Timp = TemperatureC[1]
    Tbare = TemperatureC[2]
    Tveg = TemperatureC[3]
    Ttree = TemperatureC[6]
    Tcanyon = TemperatureC[9]
    Twsun = TemperatureC[4]
    Twshade = TemperatureC[5]
    qcanyon = TemperatureC[10]
    Tveg_tm1 = TempVec_ittm.TGroundVeg
    Ttree_tm1 = TempVec_ittm.TTree
    Tcanyon_tm1 = TempVec_ittm.TCanyon

    Tveg_tm1 = TempVec_ittm.TGroundVeg;
    Ttree_tm1 = TempVec_ittm.TTree;
    Tcanyon_tm1 = TempVec_ittm.TCanyon;

    Tatm = MeteoData.Tatm;
    Pre = MeteoData.Pre;
    ea = MeteoData.ea;
    Zatm = MeteoData.Zatm;
    Uatm = MeteoData.Uatm;
    Catm_O2 = MeteoData.Catm_O2;
    Catm_CO2 = MeteoData.Catm_CO2;

    H = Gemeotry_m.Height_canyon;
    W = Gemeotry_m.Width_canyon;
    w = Gemeotry_m.wcanyon;
    Wroof = Gemeotry_m.Width_roof;
    Htree = Gemeotry_m.Height_tree;
    R_tree = Gemeotry_m.Radius_tree;
    Hcan_max = Gemeotry_m.Hcan_max;
    Hcan_std = Gemeotry_m.Hcan_std;
    wroof_norm = Gemeotry_m.wroof_norm;
    rad_tree = Gemeotry_m.radius_tree;

    fgveg = FractionsGround.fveg;
    fgbare = FractionsGround.fbare;
    fgimp = FractionsGround.fimp;

    trees = Gemeotry_m.trees;

    LAI_L = ParVegGround.LAI;
    SAI_L = ParVegGround.SAI;

    LAI_H = ParVegTree.LAI;
    SAI_H = ParVegTree.SAI;

    hc_L = ParVegGround.hc;

    hc_H = Gemeotry_m.Height_tree;

    d_leaf_L = ParVegGround.d_leaf;

    d_leaf_H = ParVegTree.d_leaf;
    Kopt_H = ParVegTree.Kopt;
    Kopt_L = ParVegGround.Kopt;
    Knit_H = ParVegTree.Knit;
    Knit_L = ParVegGround.Knit;
    Psi_sto_50_H = ParVegTree.Psi_sto_50;
    Psi_sto_50_L = ParVegGround.Psi_sto_50;
    Psi_sto_00_H = ParVegTree.Psi_sto_00;
    Psi_sto_00_L = ParVegGround.Psi_sto_00;
    CT_H = ParVegTree.CT;
    CT_L = ParVegGround.CT;
    Vmax_H = ParVegTree.Vmax;
    Vmax_L = ParVegGround.Vmax;
    DSE_H = ParVegTree.DSE;
    DSE_L = ParVegGround.DSE;
    Ha_H = ParVegTree.Ha;
    Ha_L = ParVegGround.Ha;
    FI_H = ParVegTree.FI;
    FI_L = ParVegGround.FI;
    Do_H = ParVegTree.Do;
    Do_L = ParVegGround.Do;
    a1_H = ParVegTree.a1;
    a1_L = ParVegGround.a1;
    go_H = ParVegTree.go;
    go_L = ParVegGround.go;
    e_rel_H = ParVegTree.e_rel;
    e_rel_L = ParVegGround.e_rel;
    e_relN_H = ParVegTree.e_relN;
    e_relN_L = ParVegGround.e_relN;
    gmes_H = ParVegTree.gmes;
    gmes_L = ParVegGround.gmes;
    rjv_H = ParVegTree.rjv;
    rjv_L = ParVegGround.rjv;
    mSl_H = ParVegTree.mSl;
    mSl_L = ParVegGround.mSl;
    Sl_H = ParVegTree.Sl;
    Sl_L = ParVegGround.Sl;
    SPARTREE = ParVegTree.SPARTREE;

    Pcla = ParSoilGround.Pcla;
    Psan = ParSoilGround.Psan;
    Porg = ParSoilGround.Porg;
    Kfc = ParSoilGround.Kfc;
    Phy = ParSoilGround.Phy;
    SPAR = ParSoilGround.SPAR;
    Kbot = ParSoilGround.Kbot;

    CASE_ROOT_H = ParVegTree.CASE_ROOT;
    CASE_ROOT_L = ParVegGround.CASE_ROOT;
    ZR95_H = [ParVegTree.ZR95];
    ZR95_L = [ParVegGround.ZR95];
    ZR50_H = [ParVegTree.ZR50];
    ZR50_L = [ParVegGround.ZR50];
    ZRmax_H = [ParVegTree.ZRmax];
    ZRmax_L = [ParVegGround.ZRmax];
    Zs = ParSoilGround.Zs;

    Psi_L_tm1 = SoilPotW_ittm.SoilPotWGroundVeg_L;
    Psi_H_tm1 = SoilPotW_ittm.SoilPotWGroundTot_H;

    Otm1Imp = Owater_ittm.OwGroundSoilImp;
    Otm1Bare = Owater_ittm.OwGroundSoilBare;
    Otm1Veg = Owater_ittm.OwGroundSoilVeg;
    Vtm1Imp = Vwater_ittm.VGroundSoilImp;
    Vtm1Bare = Vwater_ittm.VGroundSoilBare;
    Vtm1Veg = Vwater_ittm.VGroundSoilVeg;

    dth = ParCalculation.dth;
    row = ParCalculation.row;

    ExwatImp_tm1_H = ExWater_ittm.ExWaterGroundImp_H;
    ExwatBare_tm1_H = ExWater_ittm.ExWaterGroundBare_H;
    ExwatVeg_tm1_H = ExWater_ittm.ExWaterGroundVeg_H;
    ExwatVeg_tm1_L = ExWater_ittm.ExWaterGroundVeg_L;

    Sp_In_H = ParInterceptionTree.Sp_In;

    Sp_In_L = ParSoilGround.Sp_In;

    In_imp_tm1 = Int_ittm.IntGroundImp;
    In_bare_tm1 = Int_ittm.IntGroundBare;
    In_veg_tm1 = Int_ittm.IntGroundVegPlant;
    In_underveg_tm1 = Int_ittm.IntGroundVegGround;

    In_tree_tm1 = Int_ittm.IntTree;

    Ci_sun_H_tm1 = CiCO2Leaf_ittm.CiCO2LeafTreeSun;
    Ci_shd_H_tm1 = CiCO2Leaf_ittm.CiCO2LeafTreeShd;
    Ci_sun_L_tm1 = CiCO2Leaf_ittm.CiCO2LeafGroundVegSun;
    Ci_shd_L_tm1 = CiCO2Leaf_ittm.CiCO2LeafGroundVegShd;

    # Structural parameters
    Cimp = fgimp > 0
    Cbare = fgbare > 0
    Cveg = fgveg > 0
    Ctree = trees == 1
    Ctree025 = 4*rad_tree ≥ 0.25*w

    if isnan(rad_tree)
        rad_tree = 0.0
    end

    # Calculate thermodynamic properties
    cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    rho_atm = (Pre / (287.04 * Tatm)) * (1 - (ea / Pre) * (1 - 0.622))
    q_atm = (0.622 * ea) / (Pre - 0.378 * ea)
    L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    # Average temperatures
    Tgroundtree =
        (
            fgveg*Tveg+fgbare*Tbare+fgimp*Timp+Ctree*Ttree*(4*rad_tree)
        )/(fgveg+fgbare+fgimp+(Ctree*4*rad_tree));

    Tsurf =
        (
            fgveg*Tveg+fgbare*Tbare+fgimp*Timp+Ctree*Ttree*(4*rad_tree)+H/W*Twsun+H/W*Twshade
        )/(fgveg+fgbare+fgimp+(Ctree*4*rad_tree)+2*H/W);

    Tground = fgveg*Tveg+fgbare*Tbare+fgimp*Timp;
    Tgroundsoil = (fgveg*Tveg+fgbare*Tbare)/(fgveg+fgbare);

    # Vapor pressure calculations
    esat_T_imp = 611 * exp(17.27 * (Timp - 273.16) / (237.3 + (Timp - 273.16)))
    qsat_T_imp = (0.622 * esat_T_imp) / (Pre - 0.378 * esat_T_imp)
    esat_T_bare = 611 * exp(17.27 * (Tbare - 273.16) / (237.3 + (Tbare - 273.16)))
    qsat_T_bare = (0.622 * esat_T_bare) / (Pre - 0.378 * esat_T_bare)
    esat_T_veg = 611 * exp(17.27 * (Tveg - 273.16) / (237.3 + (Tveg - 273.16)))
    qsat_T_veg = (0.622 * esat_T_veg) / (Pre - 0.378 * esat_T_veg)
    esat_T_tree = 611 * exp(17.27 * (Ttree - 273.16) / (237.3 + (Ttree - 273.16)))
    qsat_T_tree = (0.622 * esat_T_tree) / (Pre - 0.378 * esat_T_tree)
    esat_T_canyon = 611 * exp(17.27 * (Tcanyon - 273.16) / (237.3 + (Tcanyon - 273.16)))
    qsat_T_canyon = (0.622 * esat_T_canyon) / (Pre - 0.378 * esat_T_canyon)
    e_T_canyon = qcanyon * Pre / (0.622 + 0.378 * qcanyon)
    rel_hum_canyon = e_T_canyon / esat_T_canyon
    esat_Timpbarevegtree =
        611 * exp(17.27 * (Tgroundtree - 273.16) / (237.3 + (Tgroundtree - 273.16)))

    # Parameters for stomatal resistance
    Citm1_sun_H = Ctree .* Ci_sun_H_tm1
    Citm1_shd_H = Ctree .* Ci_shd_H_tm1
    Citm1_sun_L = Cveg .* Ci_sun_L_tm1
    Citm1_shd_L = Cveg .* Ci_shd_L_tm1

    Ds_canyon = esat_T_canyon - e_T_canyon
    Oa = Catm_O2

    # Partitioning of radiation into sunlit and shaded area
    Fsun_H = Ctree * (1.0 - exp(-Kopt_H * LAI_H)) / (Kopt_H * LAI_H)
    Fsun_H = Fsun_H < 0.01 ? zero(FT) : Fsun_H
    Fsun_H = min(Fsun_H, one(FT))
    Fshd_H = Ctree * (1.0 - Fsun_H)
    PAR_sun_H = Ctree * (SWRdir_abs_tree + Fsun_H * SWRdiff_abs_tree)  # [W/m²] absorbed direct and diffuse shortwave radiation of the sunlit surface
    PAR_shd_H = Ctree * (Fshd_H * SWRdiff_abs_tree)               # [W/m²] absorbed direct and diffuse shortwave radiation of the shaded surface

    Fsun_L = Cveg * (1.0 - exp(-Kopt_L * LAI_L)) / (Kopt_L * LAI_L)
    Fsun_L = Fsun_L < 0.01 ? zero(FT) : Fsun_L
    Fsun_L = min(Fsun_L, one(FT))
    Fshd_L = Cveg*(1 - Fsun_L)
    PAR_sun_L = Cveg * (SWRdir_abs_groundveg + Fsun_L * SWRdiff_abs_groundveg)  # [W/m²] absorbed direct and diffuse shortwave radiation of the sunlit surface
    PAR_shd_L = Cveg * (Fshd_L * SWRdiff_abs_groundveg)               # [W/m²] absorbed direct and diffuse shortwave radiation of the shaded surface

    # Fraction of vegetation covered by intercepted water (Deardorff 1978)
    In_max_H = Ctree * Sp_In_H * (LAI_H + SAI_H)
    dw_H = min(1.0, (In_tree_tm1 / In_max_H)^(2/3))    # Fraction of vegetation covered by intercepted water
    In_max_L = Cveg * Sp_In_L * (LAI_L + SAI_L)
    dw_L = min(1.0, (In_veg_tm1 / In_max_L)^(2/3))    # Fraction of vegetation covered by intercepted water

    # Calculate structural parameters and wind profile in the city
    zom, zoh, zom_ground, zoh_ground, disp_h, zom_H, zom_L, zoh_H, zoh_L, d_H, d_L, zom_other = urban_roughness(
        Ctree * hc_H, Cveg * hc_L, Cbare, Cimp, false
    )

    if Ctree > 0
        dcan, zomcan, u_Hcan, u_tree, w_tree, alpha, RoughnessParameter = wind_profile_canyon(
            H,
            Htree,
            R_tree,
            W,
            Wroof,
            Kopt_H,
            LAI_H,
            Zatm,
            Uatm,
            Ctree * hc_H,
            trees,
            FT(1.5),
            zom_ground,
            Hcan_max,
            Hcan_std,
        )
    else
        u_tree = FT(0.1)
    end

    if Cveg > 0
        dcan, zomcan, u_Hcan, u_Lveg, w_Lveg, alpha, RoughnessParameter = wind_profile_canyon(
            H,
            Htree,
            R_tree,
            W,
            Wroof,
            Kopt_H,
            LAI_H,
            Zatm,
            Uatm,
            Cveg * hc_L,
            trees,
            FT(1.5),
            zom_ground,
            Hcan_max,
            Hcan_std,
        )
    else
        u_Lveg = FT(0.1)
    end

    dcan, zomcan, u_Hcan, u_Zref_und, w_Zref_und, alpha, RoughnessParameter = wind_profile_canyon(
        H,
        Htree,
        R_tree,
        W,
        Wroof,
        Kopt_H,
        LAI_H,
        Zatm,
        Uatm,
        FT(1.5),
        trees,
        FT(1.5),
        zom_ground,
        Hcan_max,
        Hcan_std,
    )

    # Calculate aerodynamic and stomata resistances
    Opt_CR = one(FT)  # [ppm] Numerical tolerance for internal CO2 computation

    rap_can, rap_Htree, rap_Htree_In, rap_2m, rap_2m_In, rap_Zp3, rap_Zp3_In, u_Hcan, u_Htree, u_Zp2, u_Zp3, uref_und, alpha = in_canyon_aerodynamic_resistance(
        Uatm,
        Zatm,
        Tcanyon-FT(273.15),
        Tsurf-FT(273.15),
        Hcan_max,
        H,
        dcan,
        zomcan,
        FT(1.5),
        zom_ground,
        hc_H,
        FT(2),
        FT(2),
        Pre,
        e_T_canyon,
        RoughnessParameter,
    )

    if Ctree == 1 && LAI_H > 0
        rb_H = leaf_boundary_resistance(
            u_tree, Ttree-FT(273.15), Tcanyon-FT(273.15), d_leaf_H, alpha
        )

        if RESPreCalc
            rs_sun_H = rsTreePreCalc.rs_sun_H
            rs_shd_H = rsTreePreCalc.rs_shd_H
            Ci_sun_H = rsTreePreCalc.Ci_sun_H
            Ci_shd_H = rsTreePreCalc.Ci_shd_H
        else
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
                Ttree_tm1-FT(273.15),
                Tcanyon_tm1-FT(273.15),
                Pre/FT(100),
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
    else
        rb_H = Inf
        rs_sun_H = Inf
        rs_shd_H = Inf
        Ci_sun_H = zero(FT)
        Ci_shd_H = zero(FT)
    end

    if Cveg == 1 && LAI_L > 0
        rb_L = leaf_boundary_resistance(
            u_Lveg, Tveg-FT(273.15), Tcanyon-FT(273.15), d_leaf_L, alpha
        )

        if RESPreCalc
            rs_sun_L = rsGroundPreCalc.rs_sun_L
            rs_shd_L = rsGroundPreCalc.rs_shd_L
            Ci_sun_L = rsGroundPreCalc.Ci_sun_L
            Ci_shd_L = rsGroundPreCalc.Ci_shd_L
        else
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
                Tveg_tm1-FT(273.15),
                Tcanyon_tm1-FT(273.15),
                Pre/FT(100),
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
    else
        rb_L = Inf
        rs_sun_L = Inf
        rs_shd_L = Inf
        Ci_sun_L = zero(FT)
        Ci_shd_L = zero(FT)
    end

    # Calculate soil parameters and resistances
    _, _, _, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, _, _, RfH_Zs, RfL_Zs, _, _, _, _, _, _, _, _, _, _ = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT_H,
        CASE_ROOT_L,
        ZR95_H,
        ZR95_L,
        ZR50_H,
        ZR50_L,
        ZRmax_H,
        ZRmax_L,
        Zs,
    )

    r_soil_bare, _, alp_soil_bare = soil_resistance(
        Tbare-FT(273.15),
        Pre/FT(100),
        u_Zref_und,
        e_T_canyon,
        In_bare_tm1,
        Otm1Bare[1],
        Ks_Zs[1],
        Osat[1],
        Ohy[1],
        L[1],
        Pe[1],
        O33[1],
        alpVG[1],
        nVG[1],
        SPAR,
    )

    r_soil_veg, _, alp_soil_veg = soil_resistance(
        Tveg-FT(273.15),
        Pre/FT(100),
        u_Zref_und,
        e_T_canyon,
        In_bare_tm1,
        Otm1Veg[1],
        Ks_Zs[1],
        Osat[1],
        Ohy[1],
        L[1],
        Pe[1],
        O33[1],
        alpVG[1],
        nVG[1],
        SPAR,
    )

    # Calculate turbulent fluxes
    Himp = Cimp * (cp_atm * rho_atm * (Timp - Tcanyon) / rap_can)
    Hbare = Cbare * (cp_atm * rho_atm * (Tbare - Tcanyon) / rap_can)
    Hveg = Cveg * (cp_atm * rho_atm * (Tveg - Tcanyon) / (rb_L/(2*(LAI_L+SAI_L)) + rap_can))
    Htree =
        Ctree *
        (cp_atm * rho_atm * (Ttree - Tcanyon) / (rb_H/(2*(LAI_H+SAI_H)) + rap_Htree_In))

    # Calculate potential evaporation rates
    Eimp_pot = Cimp * (rho_atm * (qsat_T_imp - qcanyon) / rap_can)
    Ebare_soil_pot =
        Cbare *
        (rho_atm * (alp_soil_bare * qsat_T_bare - qcanyon) / (rap_can + r_soil_bare))
    Eveg_int_pot =
        Cveg * (rho_atm * (qsat_T_veg - qcanyon) / (rb_L/((LAI_L+SAI_L)*dw_L) + rap_can))
    Eveg_soil_pot =
        Cveg * (rho_atm * (alp_soil_veg * qsat_T_veg - qcanyon) / (rap_can + r_soil_veg))

    TEveg_sun_pot =
        Cveg * (
            rho_atm * (qsat_T_veg - qcanyon) /
            (rb_L/((LAI_L)*Fsun_L*(1-dw_L)) + rap_can + rs_sun_L/((LAI_L)*Fsun_L*(1-dw_L)))
        )
    TEveg_shd_pot =
        Cveg * (
            rho_atm * (qsat_T_veg - qcanyon) /
            (rb_L/((LAI_L)*Fshd_L*(1-dw_L)) + rap_can + rs_shd_L/((LAI_L)*Fshd_L*(1-dw_L)))
        )
    TEveg_pot = Cveg * (TEveg_sun_pot + TEveg_shd_pot)

    Etree_int_pot =
        Ctree *
        (rho_atm * (qsat_T_tree - qcanyon) / (rb_H/((LAI_H+SAI_H)*dw_H) + rap_Htree_In))
    TEtree_sun_pot =
        Ctree * (
            rho_atm * (qsat_T_tree - qcanyon) / (
                rb_H/((LAI_H)*Fsun_H*(1-dw_H)) +
                rap_Htree_In +
                rs_sun_H/((LAI_H)*Fsun_H*(1-dw_H))
            )
        )
    TEtree_shd_pot =
        Ctree * (
            rho_atm * (qsat_T_tree - qcanyon) / (
                rb_H/((LAI_H)*Fshd_H*(1-dw_H)) +
                rap_Htree_In +
                rs_shd_H/((LAI_H)*Fsun_H*(1-dw_H))
            )
        )
    TEtree_pot = Ctree * (TEtree_sun_pot + TEtree_shd_pot)

    # Water limitations for interception and ponding
    Eimp = min(Eimp_pot, (In_imp_tm1/(1000*dth*3600)*row))
    Ebare_pond = min(Ebare_soil_pot, (In_bare_tm1/(1000*dth*3600)*row))
    Ebare_soil_pot = Ebare_soil_pot - Ebare_pond
    Eveg_int = min(Eveg_int_pot, (In_veg_tm1/(1000*dth*3600)*row))
    Eveg_pond = min(Eveg_soil_pot, (In_underveg_tm1/(1000*dth*3600)*row))
    Eveg_soil_pot = Eveg_soil_pot - Eveg_pond
    Etree_int = min(Etree_int_pot, (In_tree_tm1/(1000*dth*3600)*row))

    # Water limitations for soil evaporation
    VavailImp_tm1 = (Vtm1Imp ./ dth) .* (row/3600/1000)
    VavailBare_tm1 = (Vtm1Bare ./ dth) .* (row/3600/1000)
    VavailVeg_tm1 = (Vtm1Veg ./ dth) .* (row/3600/1000)
    Ebare_soil = min(Ebare_soil_pot, VavailBare_tm1[1])
    Eveg_soil = min(Eveg_soil_pot, VavailVeg_tm1[1])
    VavailBare_tm1[1] -= Ebare_soil
    VavailVeg_tm1[1] -= Eveg_soil

    # Water limitations for vegetation transpiration
    RfH_Zs_Imp = Vector{FT}(undef, length(Zs)-1)
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, RfH_Zs_ImpL2, _, _, _, _, _, _, _, _, _, _, _ = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT_H,
        CASE_ROOT_L,
        ZR95_H,
        ZR95_L,
        ZR50_H,
        ZR50_L,
        ZRmax_H,
        ZRmax_L,
        Zs[3:end],
    )
    RfH_Zs_Imp[1:2] .= zero(FT)
    RfH_Zs_Imp[3:end] = RfH_Zs_ImpL2

    # Calculate crown areas and available water
    if SPARTREE == 1 # Tree roots can access all water in soil
        Ccrown = [Cveg*fgveg, Cbare*fgbare, Cimp*fgimp, Ctree*4*rad_tree]

        Vavail_Veg_tm1_L = Ccrown[1]/(Ccrown[1]+Ccrown[1]*Ccrown[4]) .* VavailVeg_tm1
        Vavail_Veg_tm1_H =
            Ccrown[1]*Ccrown[4]/(Ccrown[1]+Ccrown[1]*Ccrown[4]) .* VavailVeg_tm1
        Vavail_Bare_tm1_H = Ccrown[2]*Ccrown[4]/(Ccrown[2]*Ccrown[4]) .* VavailBare_tm1
        Vavail_Imp_tm1_H = Ccrown[3]*Ccrown[4]/(Ccrown[3]*Ccrown[4]) .* VavailImp_tm1

    else # SPARTREE == 2: Limited tree root access
        if (4*rad_tree) ≤ (fgveg+fgbare)
            Ccrown = [Cveg*fgveg, Cbare*fgbare, Cimp*fgimp, Ctree*4*rad_tree/(fgveg+fgbare)]

            Vavail_Veg_tm1_L = Ccrown[1]/(Ccrown[1]+Ccrown[1]*Ccrown[4]) .* VavailVeg_tm1
            Vavail_Veg_tm1_H =
                Ccrown[1]*Ccrown[4]/(Ccrown[1]+Ccrown[1]*Ccrown[4]) .* VavailVeg_tm1
            Vavail_Bare_tm1_H = Ccrown[2]*Ccrown[4]/(Ccrown[2]*Ccrown[4]) .* VavailBare_tm1
            Vavail_Imp_tm1_H = zero(VavailImp_tm1)

        else # (4*rad_tree) > (fgveg+fgbare)
            Ccrown = [
                Cveg*fgveg,
                Cbare*fgbare,
                Cimp*fgimp,
                Ctree*one(FT),
                Ctree*((4*rad_tree)-(fgveg+fgbare))/fgimp,
            ]

            Vavail_Veg_tm1_L = Ccrown[1]/(Ccrown[1]+Ccrown[1]*Ccrown[4]) .* VavailVeg_tm1
            Vavail_Veg_tm1_H =
                Ccrown[1]*Ccrown[4]/(Ccrown[1]+Ccrown[1]*Ccrown[4]) .* VavailVeg_tm1
            Vavail_Bare_tm1_H = Ccrown[2]*Ccrown[4]/(Ccrown[2]*Ccrown[4]) .* VavailBare_tm1
            Vavail_Imp_tm1_H = Ccrown[3]*Ccrown[5]/(Ccrown[3]*Ccrown[5]) .* VavailImp_tm1
        end
    end

    # Calculate minimum available and extractable water
    Vavail_Veg_tm1_L = sum(
        NaNMath.min.(Vavail_Veg_tm1_L .* RfL_Zs, ExwatVeg_tm1_L*(row/3600/1000))
    )
    Vavail_Veg_tm1_H = sum(
        NaNMath.min.(Vavail_Veg_tm1_H .* RfH_Zs, ExwatVeg_tm1_H*(row/3600/1000))
    )
    Vavail_Bare_tm1_H = sum(
        NaNMath.min.(Vavail_Bare_tm1_H .* RfH_Zs, ExwatBare_tm1_H*(row/3600/1000))
    )
    Vavail_Imp_tm1_H = sum(
        NaNMath.min.(Vavail_Imp_tm1_H .* RfH_Zs_Imp, ExwatImp_tm1_H*(row/3600/1000))
    )

    # Water limitation for transpiration
    Vavail_plant_tm1_H =
        (
            Ccrown[1]*Vavail_Veg_tm1_H +
            Ccrown[2]*Vavail_Bare_tm1_H +
            Ccrown[3]*Vavail_Imp_tm1_H
        )/Ccrown[4]
    Vavail_plant_tm1_L = Vavail_Veg_tm1_L
    TEveg = min(TEveg_pot, Vavail_plant_tm1_L)
    TEtree = min(TEtree_pot, Vavail_plant_tm1_H)

    # Calculate final evapotranspiration rates
    Ebare = Cbare * (Ebare_pond + Ebare_soil)
    Eveg = Cveg * (Eveg_int + Eveg_pond + Eveg_soil + TEveg)
    Etree = Ctree * (Etree_int + TEtree)

    # Calculate latent heat fluxes
    LEimp = Cimp * (L_heat * Eimp)
    LEbare_pond = Cbare * (L_heat * Ebare_pond)
    LEbare_soil = Cbare * (L_heat * Ebare_soil)
    LEbare = Cbare * (LEbare_pond + LEbare_soil)
    LEveg_int = Cveg * (L_heat * Eveg_int)
    LEveg_pond = Cveg * (L_heat * Eveg_pond)
    LEveg_soil = Cveg * (L_heat * Eveg_soil)
    LTEveg = Cveg * (L_heat * TEveg)
    LEveg = Cveg * (LEveg_int + LEveg_pond + LEveg_soil + LTEveg)
    LEtree_int = Ctree * (L_heat * Etree_int)
    LTEtree = Ctree * (L_heat * TEtree)
    LEtree = Ctree * (LEtree_int + LTEtree)

    return (;
        Himp,
        Hbare,
        Hveg,
        Htree,
        Eimp,
        Ebare_pond,
        Ebare_soil,
        Eveg_int,
        Eveg_pond,
        Eveg_soil,
        TEveg,
        Etree_int,
        TEtree,
        Ebare,
        Eveg,
        Etree,
        LEimp,
        LEbare_pond,
        LEbare_soil,
        LEveg_int,
        LEveg_pond,
        LEveg_soil,
        LTEveg,
        LEtree_int,
        LTEtree,
        LEbare,
        LEveg,
        LEtree,
        Ci_sun_H,
        Ci_shd_H,
        Ci_sun_L,
        Ci_shd_L,
        rap_can,
        rap_Htree_In,
        rb_H,
        rb_L,
        r_soil_bare,
        r_soil_veg,
        alp_soil_bare,
        alp_soil_veg,
        rs_sun_L,
        rs_shd_L,
        rs_sun_H,
        rs_shd_H,
        u_Hcan,
        u_Zref_und,
        Fsun_L,
        Fshd_L,
        dw_L,
    )
end
