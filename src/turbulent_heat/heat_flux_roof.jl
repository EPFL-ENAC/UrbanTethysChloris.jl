"""
    heat_flux_roof(
        TemperatureR::AbstractVector{FT},
        TempVec_ittm::NamedTuple,
        MeteoData::NamedTuple,
        HumidityAtm::NamedTuple,
        ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
        Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
        ParCalculation::NamedTuple,
        SoilPotW_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        Vwater_ittm::NamedTuple,
        ExWater_ittm::NamedTuple,
        Int_ittm::NamedTuple,
        CiCO2Leaf_ittm::NamedTuple,
        SWRabs_dir::FT,
        SWRabs_diff::FT,
        RESPreCalc::Bool,
        rsRoofPreCalc::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate sensible and latent heat fluxes for roof surfaces.

# Arguments
- `TemperatureR`: Temperature vector [K]
- `TempVec_ittm`: Temperature variables from previous timestep
- `MeteoData`: Meteorological data
- `HumidityAtm`: Atmospheric humidity parameters
- `ParVegRoof`: Roof vegetation parameters
- `FractionsRoof`: Roof surface fractions
- `Gemeotry_m`: Urban geometry parameters
- `ParSoilRoof`: Roof soil parameters
- `ParCalculation`: Calculation parameters
- `SoilPotW_ittm`: Soil water potential from previous timestep
- `Owater_ittm`: Soil water content from previous timestep
- `Vwater_ittm`: Soil water volume from previous timestep
- `ExWater_ittm`: Extractable water from previous timestep
- `Int_ittm`: Interception from previous timestep
- `CiCO2Leaf_ittm`: Leaf CO2 concentration from previous timestep
- `SWRabs_dir`: Direct shortwave radiation [W/m²]
- `SWRabs_diff`: Diffuse shortwave radiation [W/m²]
- `RESPreCalc`: Use pre-calculated resistances
- `rsRoofPreCalc`: Pre-calculated resistance parameters

# Returns
- `Hroof_imp::FT`: Sensible heat flux from impervious roof [W/m²]
- `Hroof_veg::FT`: Sensible heat flux from vegetated roof [W/m²]
- `Eroof_imp::FT`: Water vapor flux from impervious roof [kg/m²s]
- `Eroof_veg::FT`: Water vapor flux from vegetated roof [kg/m²s]
- `Eroof_ground::FT`: Ground water vapor flux from vegetated roof [kg/m²s]
- `Eroof_soil::FT`: Soil water vapor flux from vegetated roof [kg/m²s]
- `TEroof_veg::FT`: Transpiration flux from vegetated roof [kg/m²s]
- `LEroof_imp::FT`: Latent heat flux from impervious roof [W/m²]
- `LEroof_veg::FT`: Latent heat flux from vegetated roof [W/m²]
- `LEroof_ground::FT`: Ground latent heat flux from vegetated roof [W/m²]
- `LEroof_soil::FT`: Soil latent heat flux from vegetated roof [W/m²]
- `LTEroof_veg::FT`: Latent heat of transpiration from vegetated roof [W/m²]
- `Ci_sun_roof::FT`: Sunlit leaf internal CO2 concentration [μmol/mol]
- `Ci_shd_roof::FT`: Shaded leaf internal CO2 concentration [μmol/mol]
- `ra::FT`: Atmospheric resistance [s/m]
- `rb::FT`: Boundary layer resistance [s/m]
- `rap_L::FT`: Roof aerodynamic resistance [s/m]
- `r_soil::FT`: Soil resistance [s/m]
- `rs_sun::FT`: Sunlit stomatal resistance [s/m]
- `rs_shd::FT`: Shaded stomatal resistance [s/m]
"""
function heat_flux_roof(
    TemperatureR::AbstractVector{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    HumidityAtm::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParCalculation::NamedTuple,
    SoilPotW_ittm::ModelComponents.ModelVariables.SoilPotW{FT},
    Owater_ittm::ModelComponents.ModelVariables.Owater{FT,MR,MG},
    Vwater_ittm::ModelComponents.ModelVariables.Vwater{FT,MR,MG},
    ExWater_ittm::ModelComponents.ModelVariables.ExWater{FT,MR,MG},
    Int_ittm::ModelComponents.ModelVariables.Interception{FT},
    CiCO2Leaf_ittm::ModelComponents.ModelVariables.CiCO2Leaf{FT},
    SWRabs_dir::FT,
    SWRabs_diff::FT,
    RESPreCalc::Bool,
    rsRoofPreCalc::NamedTuple,
) where {FT<:AbstractFloat,MR,MG}
    # Extract temperatures
    Troof_imp = TemperatureR[1]
    Troof_veg = TemperatureR[2]
    Troof_veg_tm1 = TempVec_ittm.TRoofVeg

    # Extract meteorological data
    Tatm = MeteoData.Tatm
    Pre = MeteoData.Pre
    ea = MeteoData.ea
    Zatm = MeteoData.Zatm
    Uatm = MeteoData.Uatm
    Catm_O2 = MeteoData.Catm_O2
    Catm_CO2 = MeteoData.Catm_CO2
    esat_Tatm = HumidityAtm.AtmVapourPreSat

    # Extract geometry parameters
    H = Gemeotry_m.Height_canyon

    # Extract vegetation parameters
    fveg_roof = FractionsRoof.fveg
    fimp_roof = FractionsRoof.fimp
    LAI_roof = ParVegRoof.LAI
    SAI_roof = ParVegRoof.SAI
    hc_roof = ParVegRoof.hc
    d_leaf_roof = ParVegRoof.d_leaf
    Kopt_roof = ParVegRoof.Kopt
    Knit_roof = ParVegRoof.Knit
    Psi_sto_50_roof = ParVegRoof.Psi_sto_50
    Psi_sto_00_roof = ParVegRoof.Psi_sto_00
    CT_roof = ParVegRoof.CT
    Vmax_roof = ParVegRoof.Vmax
    DSE_roof = ParVegRoof.DSE
    Ha_roof = ParVegRoof.Ha
    FI_roof = ParVegRoof.FI
    Do_roof = ParVegRoof.Do
    a1_roof = ParVegRoof.a1
    go_roof = ParVegRoof.go
    e_rel_roof = ParVegRoof.e_rel
    e_relN_roof = ParVegRoof.e_relN
    gmes_roof = ParVegRoof.gmes
    rjv_roof = ParVegRoof.rjv
    mSl_roof = ParVegRoof.mSl
    Sl_roof = ParVegRoof.Sl

    # Extract root parameters
    CASE_ROOT = ParVegRoof.CASE_ROOT
    ZR95 = ParVegRoof.ZR95
    ZR50 = ParVegRoof.ZR50
    ZRmax = ParVegRoof.ZRmax

    # Extract soil parameters
    Pcla = ParSoilRoof.Pcla
    Psan = ParSoilRoof.Psan
    Porg = ParSoilRoof.Porg
    Kfc = ParSoilRoof.Kfc
    Phy = ParSoilRoof.Phy
    SPAR = ParSoilRoof.SPAR
    Kbot = ParSoilRoof.Kbot
    Sp_In_roof = ParSoilRoof.Sp_In
    Zs = ParSoilRoof.Zs

    # Extract calculation parameters
    dth = ParCalculation.dth
    row = ParCalculation.row

    # Extract previous timestep variables
    Psi_ltm1 = SoilPotW_ittm.SoilPotWRoofVeg_L
    Otm1 = Owater_ittm.OwRoofSoilVeg
    Vtm1 = Vwater_ittm.VRoofSoilVeg
    Exwat_tm1 = ExWater_ittm.ExWaterRoofVeg_L
    In_imp_tm1 = Int_ittm.IntRoofImp
    In_veg_tm1 = Int_ittm.IntRoofVegPlant
    In_ground_tm1 = Int_ittm.IntRoofVegGround
    Ci_sun_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegSun
    Ci_shd_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegShd

    ## Parameters
    # Calculate thermodynamic properties
    cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    rho_atm = (Pre / (287.04 * Tatm)) * (1 - (ea / Pre) * (1 - 0.622))
    q_atm = 0.622 * ea / (Pre - 0.378 * ea)
    L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    # Vapor pressure calculations
    esat_T_rimp = 611 * exp(17.27 * (Troof_imp - 273.16) / (237.3 + (Troof_imp - 273.16)))
    qsat_T_rimp = (0.622 * esat_T_rimp) / (Pre - 0.378 * esat_T_rimp)
    esat_T_rveg = 611 * exp(17.27 * (Troof_veg - 273.16) / (237.3 + (Troof_veg - 273.16)))
    qsat_T_rveg = (0.622 * esat_T_rveg) / (Pre - 0.378 * esat_T_rveg)
    Troof = fveg_roof * Troof_veg + fimp_roof * Troof_imp
    esat_roof = fveg_roof * esat_T_rveg + fimp_roof * esat_T_rimp

    # Radiation partitioning
    Fsun = (1.0 - exp(-Kopt_roof * LAI_roof)) / (Kopt_roof * LAI_roof)
    Fsun = clamp(Fsun, 0.0, 1.0)
    Fshd = 1 - Fsun
    PAR_sun = SWRabs_dir + Fsun * SWRabs_diff
    PAR_shd = Fshd * SWRabs_diff

    # Interception calculations
    In_max_veg = Sp_In_roof * (LAI_roof + SAI_roof)
    dw_veg_roof = min(1.0, (In_veg_tm1 / In_max_veg)^(2/3))

    ## Calculation of resistances
    # Calculate roughness parameters
    if fveg_roof == 0
        hc_roof = zero(FT)
        LAI_roof = zero(FT)
        d_leaf_roof = zero(FT)
    end
    Croof = fimp_roof == 0 ? false : true

    zom, zoh, _, _, disp_h = urban_roughness(hc_roof, zero(FT), false, false, Croof)

    # Wind profile and resistances
    u_ref_und, u_Hveg = wind_profile_roof(H, Zatm, Uatm, zom, disp_h, hc_roof, FT(2))

    ra = aerodynamic_resistance(
        Tatm - 273.15,
        Troof - 273.15,
        Pre/100,
        Zatm - H,
        disp_h,
        zom,
        zoh,
        Uatm,
        ea,
        esat_roof,
    )

    rb = leaf_boundary_resistance(
        u_Hveg,
        Troof_veg - 273.15,
        Tatm - 273,
        hc_roof,
        d_leaf_roof,
        LAI_roof,
        Zatm,
        disp_h,
        zom,
    )

    if RESPreCalc
        rs_sun = rsRoofPreCalc.rs_sun
        rs_shd = rsRoofPreCalc.rs_shd
        Ci_sun_roof = rsRoofPreCalc.Ci_sun
        Ci_shd_roof = rsRoofPreCalc.Ci_shd
    else
        # Parameters for stomata resistance
        Citm1_sun = Ci_sun_tm1      # Leaf Interior CO2 concentration [umolCO2/mol]
        Citm1_shd = Ci_shd_tm1      # Leaf Interior CO2 concentration [umolCO2/mol]
        Ds_atm = esat_Tatm - ea     # Vapor Pressure Deficit [Pa]
        Oa = Catm_O2                # Intercellular Partial Pressure Oxygen [umolO2/mol]

        # Stomatal resistance calculation
        if LAI_roof > 0
            rs_sun, rs_shd, Ci_sun_roof, Ci_shd_roof, _, _, _, _, _ = canopy_resistance_an_evolution(
                PAR_sun,
                PAR_shd,
                LAI_roof,
                Kopt_roof,
                Knit_roof,
                Fsun,
                Fshd,
                Citm1_sun,
                Citm1_shd,
                Catm_CO2,
                ra,
                rb,
                Troof_veg_tm1-273.15,
                Tatm-273.15,
                Pre/100,
                Ds_atm,
                Psi_ltm1,
                Psi_sto_50_roof,
                Psi_sto_00_roof,
                CT_roof,
                Vmax_roof,
                DSE_roof,
                Ha_roof,
                FI_roof,
                Oa,
                Do_roof,
                a1_roof,
                go_roof,
                e_rel_roof,
                e_relN_roof,
                gmes_roof,
                rjv_roof,
                mSl_roof,
                Sl_roof,
                FT(1),  # Numerical tolerance
            )
        else
            rs_sun = Inf
            rs_shd = Inf
            Ci_sun_roof = zero(FT)
            Ci_shd_roof = zero(FT)
        end
    end

    # Soil parameters and resistance
    _, _, _, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, _, _, _, Rf_Zs, _, _ = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT,
        CASE_ROOT,
        ZR95,
        ZR95,
        ZR50,
        ZR50,
        ZRmax,
        ZRmax,
        Zs,
    )

    r_soil, _, alp_soil = soil_resistance(
        Troof_veg-273.15,
        Pre/100,
        u_ref_und,
        ea,
        In_ground_tm1,
        Otm1[1],
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

    ## Calculatation of turbulent fluxes
    Hroof_imp = cp_atm * rho_atm * (Troof_imp - Tatm) / ra
    Hroof_veg = cp_atm * rho_atm * (Troof_veg - Tatm) / (rb/(2*(LAI_roof + SAI_roof)) + ra)

    Eroof_imp_pot = rho_atm * (qsat_T_rimp - q_atm) / ra
    Eroof_veg_pot =
        rho_atm * (qsat_T_rveg - q_atm) / (rb/((LAI_roof + SAI_roof) * dw_veg_roof) + ra)
    Eroof_soil_pot = rho_atm * alp_soil * (qsat_T_rveg - q_atm) / (ra + r_soil)

    TEroof_veg_sun_pot =
        rho_atm * (qsat_T_rveg - q_atm) / (
            rb/((LAI_roof) * (1-dw_veg_roof)) +
            ra +
            rs_sun/((LAI_roof) * Fsun * (1-dw_veg_roof))
        )
    TEroof_veg_shd_pot =
        rho_atm * (qsat_T_rveg - q_atm) / (
            rb/((LAI_roof) * (1-dw_veg_roof)) +
            ra +
            rs_shd/((LAI_roof) * Fshd * (1-dw_veg_roof))
        )
    TEroof_veg_pot = TEroof_veg_sun_pot + TEroof_veg_shd_pot

    # Ensure evapotranspiration doesn't exceed available water
    Eroof_imp = min(Eroof_imp_pot, (In_imp_tm1/(1000*dth*3600)*row))
    Eroof_veg = min(Eroof_veg_pot, (In_veg_tm1/(1000*dth*3600)*row))
    Eroof_ground = min(Eroof_soil_pot, (In_ground_tm1/(1000*dth*3600)*row))
    Eroof_soil_pot = Eroof_soil_pot - Eroof_ground

    # Water available for Transpiration and Evaporation
    Vavail_tm1 = (Vtm1 ./ dth) .* (row/3600/1000)
    Eroof_soil = min(Eroof_soil_pot, Vavail_tm1[1])
    Vavail_tm1[1] = Vavail_tm1[1] - Eroof_soil
    Vavail_plant_tm1 = min(Vavail_tm1 .* vec(Rf_Zs), Exwat_tm1 .* (row/3600/1000))
    Vavail_plant_tm1 = sum(Vavail_plant_tm1)
    TEroof_veg = min(TEroof_veg_pot, Vavail_plant_tm1)

    # Apply vegetation fraction conditions
    if fveg_roof == 0
        Hroof_veg = zero(FT)
        Eroof_veg = zero(FT)
        Eroof_ground = zero(FT)
        Eroof_soil = zero(FT)
        TEroof_veg = zero(FT)
    end

    if fimp_roof == 0
        Hroof_imp = zero(FT)
        Eroof_imp = zero(FT)
    end

    # Calculate latent heat fluxes
    LEroof_imp = L_heat * Eroof_imp
    LEroof_veg = L_heat * Eroof_veg
    LEroof_ground = L_heat * Eroof_ground
    LEroof_soil = L_heat * Eroof_soil
    LTEroof_veg = L_heat * TEroof_veg

    rap_L = FT(NaN)

    return Hroof_imp,
    Hroof_veg,
    Eroof_imp,
    Eroof_veg,
    Eroof_ground,
    Eroof_soil,
    TEroof_veg,
    LEroof_imp,
    LEroof_veg,
    LEroof_ground,
    LEroof_soil,
    LTEroof_veg,
    Ci_sun_roof,
    Ci_shd_roof,
    ra,
    rb,
    rap_L,
    r_soil,
    rs_sun,
    rs_shd
end
function heat_flux_roof(
    TemperatureR::AbstractVector{FT},
    TempVec_ittm::NamedTuple,
    MeteoData::NamedTuple,
    HumidityAtm::NamedTuple,
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    FractionsRoof::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    ParSoilRoof::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParCalculation::NamedTuple,
    SoilPotW_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    Vwater_ittm::NamedTuple,
    ExWater_ittm::NamedTuple,
    Int_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    SWRabs_dir::FT,
    SWRabs_diff::FT,
    RESPreCalc::Bool,
    rsRoofPreCalc::NamedTuple,
) where {FT<:AbstractFloat}
    # Extract temperatures
    Troof_imp = TemperatureR[1]
    Troof_veg = TemperatureR[2]
    Troof_veg_tm1 = TempVec_ittm.TRoofVeg

    # Extract meteorological data
    Tatm = MeteoData.Tatm
    Pre = MeteoData.Pre
    ea = MeteoData.ea
    Zatm = MeteoData.Zatm
    Uatm = MeteoData.Uatm
    Catm_O2 = MeteoData.Catm_O2
    Catm_CO2 = MeteoData.Catm_CO2
    esat_Tatm = HumidityAtm.AtmVapourPreSat

    # Extract geometry parameters
    H = Gemeotry_m.Height_canyon

    # Extract vegetation parameters
    fveg_roof = FractionsRoof.fveg
    fimp_roof = FractionsRoof.fimp
    LAI_roof = ParVegRoof.LAI
    SAI_roof = ParVegRoof.SAI
    hc_roof = ParVegRoof.hc
    d_leaf_roof = ParVegRoof.d_leaf
    Kopt_roof = ParVegRoof.Kopt
    Knit_roof = ParVegRoof.Knit
    Psi_sto_50_roof = ParVegRoof.Psi_sto_50
    Psi_sto_00_roof = ParVegRoof.Psi_sto_00
    CT_roof = ParVegRoof.CT
    Vmax_roof = ParVegRoof.Vmax
    DSE_roof = ParVegRoof.DSE
    Ha_roof = ParVegRoof.Ha
    FI_roof = ParVegRoof.FI
    Do_roof = ParVegRoof.Do
    a1_roof = ParVegRoof.a1
    go_roof = ParVegRoof.go
    e_rel_roof = ParVegRoof.e_rel
    e_relN_roof = ParVegRoof.e_relN
    gmes_roof = ParVegRoof.gmes
    rjv_roof = ParVegRoof.rjv
    mSl_roof = ParVegRoof.mSl
    Sl_roof = ParVegRoof.Sl

    # Extract root parameters
    CASE_ROOT = ParVegRoof.CASE_ROOT
    ZR95 = ParVegRoof.ZR95
    ZR50 = ParVegRoof.ZR50
    ZRmax = ParVegRoof.ZRmax

    # Extract soil parameters
    Pcla = ParSoilRoof.Pcla
    Psan = ParSoilRoof.Psan
    Porg = ParSoilRoof.Porg
    Kfc = ParSoilRoof.Kfc
    Phy = ParSoilRoof.Phy
    SPAR = ParSoilRoof.SPAR
    Kbot = ParSoilRoof.Kbot
    Sp_In_roof = ParSoilRoof.Sp_In
    Zs = ParSoilRoof.Zs

    # Extract calculation parameters
    dth = ParCalculation.dth
    row = ParCalculation.row

    # Extract previous timestep variables
    Psi_ltm1 = SoilPotW_ittm.SoilPotWRoofVeg_L
    Otm1 = Owater_ittm.OwRoofSoilVeg
    Vtm1 = Vwater_ittm.VRoofSoilVeg
    Exwat_tm1 = ExWater_ittm.ExWaterRoofVeg_L
    In_imp_tm1 = Int_ittm.IntRoofImp
    In_veg_tm1 = Int_ittm.IntRoofVegPlant
    In_ground_tm1 = Int_ittm.IntRoofVegGround
    Ci_sun_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegSun
    Ci_shd_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegShd

    ## Parameters
    # Calculate thermodynamic properties
    cp_atm = 1005 + ((Tatm - 273.15) + 23.15)^2 / 3364
    rho_atm = (Pre / (287.04 * Tatm)) * (1 - (ea / Pre) * (1 - 0.622))
    q_atm = 0.622 * ea / (Pre - 0.378 * ea)
    L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    # Vapor pressure calculations
    esat_T_rimp = 611 * exp(17.27 * (Troof_imp - 273.16) / (237.3 + (Troof_imp - 273.16)))
    qsat_T_rimp = (0.622 * esat_T_rimp) / (Pre - 0.378 * esat_T_rimp)
    esat_T_rveg = 611 * exp(17.27 * (Troof_veg - 273.16) / (237.3 + (Troof_veg - 273.16)))
    qsat_T_rveg = (0.622 * esat_T_rveg) / (Pre - 0.378 * esat_T_rveg)
    Troof = fveg_roof * Troof_veg + fimp_roof * Troof_imp
    esat_roof = fveg_roof * esat_T_rveg + fimp_roof * esat_T_rimp

    # Radiation partitioning
    Fsun = (1.0 - exp(-Kopt_roof * LAI_roof)) / (Kopt_roof * LAI_roof)
    Fsun = clamp(Fsun, 0.0, 1.0)
    Fshd = 1 - Fsun
    PAR_sun = SWRabs_dir + Fsun * SWRabs_diff
    PAR_shd = Fshd * SWRabs_diff

    # Interception calculations
    In_max_veg = Sp_In_roof * (LAI_roof + SAI_roof)
    dw_veg_roof = min(1.0, (In_veg_tm1 / In_max_veg)^(2/3))

    ## Calculation of resistances
    # Calculate roughness parameters
    if fveg_roof == 0
        hc_roof = zero(FT)
        LAI_roof = zero(FT)
        d_leaf_roof = zero(FT)
    end
    Croof = fimp_roof == 0 ? false : true

    zom, zoh, _, _, disp_h = urban_roughness(hc_roof, zero(FT), false, false, Croof)

    # Wind profile and resistances
    u_ref_und, u_Hveg = wind_profile_roof(H, Zatm, Uatm, zom, disp_h, hc_roof, FT(2))

    ra = aerodynamic_resistance(
        Tatm - 273.15,
        Troof - 273.15,
        Pre/100,
        Zatm - H,
        disp_h,
        zom,
        zoh,
        Uatm,
        ea,
        esat_roof,
    )

    rb = leaf_boundary_resistance(
        u_Hveg,
        Troof_veg - 273.15,
        Tatm - 273,
        hc_roof,
        d_leaf_roof,
        LAI_roof,
        Zatm,
        disp_h,
        zom,
    )

    if RESPreCalc
        rs_sun = rsRoofPreCalc.rs_sun
        rs_shd = rsRoofPreCalc.rs_shd
        Ci_sun_roof = rsRoofPreCalc.Ci_sun
        Ci_shd_roof = rsRoofPreCalc.Ci_shd
    else
        # Parameters for stomata resistance
        Citm1_sun = Ci_sun_tm1      # Leaf Interior CO2 concentration [umolCO2/mol]
        Citm1_shd = Ci_shd_tm1      # Leaf Interior CO2 concentration [umolCO2/mol]
        Ds_atm = esat_Tatm - ea     # Vapor Pressure Deficit [Pa]
        Oa = Catm_O2                # Intercellular Partial Pressure Oxygen [umolO2/mol]

        # Stomatal resistance calculation
        if LAI_roof > 0
            rs_sun, rs_shd, Ci_sun_roof, Ci_shd_roof, _, _, _, _, _ = canopy_resistance_an_evolution(
                PAR_sun,
                PAR_shd,
                LAI_roof,
                Kopt_roof,
                Knit_roof,
                Fsun,
                Fshd,
                Citm1_sun,
                Citm1_shd,
                Catm_CO2,
                ra,
                rb,
                Troof_veg_tm1-273.15,
                Tatm-273.15,
                Pre/100,
                Ds_atm,
                Psi_ltm1,
                Psi_sto_50_roof,
                Psi_sto_00_roof,
                CT_roof,
                Vmax_roof,
                DSE_roof,
                Ha_roof,
                FI_roof,
                Oa,
                Do_roof,
                a1_roof,
                go_roof,
                e_rel_roof,
                e_relN_roof,
                gmes_roof,
                rjv_roof,
                mSl_roof,
                Sl_roof,
                FT(1),  # Numerical tolerance
            )
        else
            rs_sun = Inf
            rs_shd = Inf
            Ci_sun_roof = zero(FT)
            Ci_shd_roof = zero(FT)
        end
    end

    # Soil parameters and resistance
    _, _, _, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, _, _, _, Rf_Zs, _, _ = soil_parameters_total(
        Pcla,
        Psan,
        Porg,
        Kfc,
        Phy,
        SPAR,
        Kbot,
        CASE_ROOT,
        CASE_ROOT,
        ZR95,
        ZR95,
        ZR50,
        ZR50,
        ZRmax,
        ZRmax,
        Zs,
    )

    r_soil, _, alp_soil = soil_resistance(
        Troof_veg-273.15,
        Pre/100,
        u_ref_und,
        ea,
        In_ground_tm1,
        Otm1[1],
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

    ## Calculatation of turbulent fluxes
    Hroof_imp = cp_atm * rho_atm * (Troof_imp - Tatm) / ra
    Hroof_veg = cp_atm * rho_atm * (Troof_veg - Tatm) / (rb/(2*(LAI_roof + SAI_roof)) + ra)

    Eroof_imp_pot = rho_atm * (qsat_T_rimp - q_atm) / ra
    Eroof_veg_pot =
        rho_atm * (qsat_T_rveg - q_atm) / (rb/((LAI_roof + SAI_roof) * dw_veg_roof) + ra)
    Eroof_soil_pot = rho_atm * alp_soil * (qsat_T_rveg - q_atm) / (ra + r_soil)

    TEroof_veg_sun_pot =
        rho_atm * (qsat_T_rveg - q_atm) / (
            rb/((LAI_roof) * (1-dw_veg_roof)) +
            ra +
            rs_sun/((LAI_roof) * Fsun * (1-dw_veg_roof))
        )
    TEroof_veg_shd_pot =
        rho_atm * (qsat_T_rveg - q_atm) / (
            rb/((LAI_roof) * (1-dw_veg_roof)) +
            ra +
            rs_shd/((LAI_roof) * Fshd * (1-dw_veg_roof))
        )
    TEroof_veg_pot = TEroof_veg_sun_pot + TEroof_veg_shd_pot

    # Ensure evapotranspiration doesn't exceed available water
    Eroof_imp = min(Eroof_imp_pot, (In_imp_tm1/(1000*dth*3600)*row))
    Eroof_veg = min(Eroof_veg_pot, (In_veg_tm1/(1000*dth*3600)*row))
    Eroof_ground = min(Eroof_soil_pot, (In_ground_tm1/(1000*dth*3600)*row))
    Eroof_soil_pot = Eroof_soil_pot - Eroof_ground

    # Water available for Transpiration and Evaporation
    Vavail_tm1 = (Vtm1 ./ dth) .* (row/3600/1000)
    Eroof_soil = min(Eroof_soil_pot, Vavail_tm1[1])
    Vavail_tm1[1] = Vavail_tm1[1] - Eroof_soil
    Vavail_plant_tm1 = min(Vavail_tm1 .* vec(Rf_Zs), Exwat_tm1 .* (row/3600/1000))
    Vavail_plant_tm1 = sum(Vavail_plant_tm1)
    TEroof_veg = min(TEroof_veg_pot, Vavail_plant_tm1)

    # Apply vegetation fraction conditions
    if fveg_roof == 0
        Hroof_veg = zero(FT)
        Eroof_veg = zero(FT)
        Eroof_ground = zero(FT)
        Eroof_soil = zero(FT)
        TEroof_veg = zero(FT)
    end

    if fimp_roof == 0
        Hroof_imp = zero(FT)
        Eroof_imp = zero(FT)
    end

    # Calculate latent heat fluxes
    LEroof_imp = L_heat * Eroof_imp
    LEroof_veg = L_heat * Eroof_veg
    LEroof_ground = L_heat * Eroof_ground
    LEroof_soil = L_heat * Eroof_soil
    LTEroof_veg = L_heat * TEroof_veg

    rap_L = FT(NaN)

    return Hroof_imp,
    Hroof_veg,
    Eroof_imp,
    Eroof_veg,
    Eroof_ground,
    Eroof_soil,
    TEroof_veg,
    LEroof_imp,
    LEroof_veg,
    LEroof_ground,
    LEroof_soil,
    LTEroof_veg,
    Ci_sun_roof,
    Ci_shd_roof,
    ra,
    rb,
    rap_L,
    r_soil,
    rs_sun,
    rs_shd
end
