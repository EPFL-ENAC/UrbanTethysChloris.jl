"""
    precalculate_stomatal_resistance_roof(
        TempVec_ittm::NamedTuple,
        MeteoData::NamedTuple,
        HumidityAtm::NamedTuple,
        ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
        SoilPotW_ittm::NamedTuple,
        CiCO2Leaf_ittm::NamedTuple,
        PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
        ra::FT,
        rb::FT
    ) where {FT<:AbstractFloat}

Calculate stomatal resistance for sunlit and shaded portions of the roof vegetation.

# Arguments
- `TempVec_ittm`: Temperature vector containing roof vegetation temperature
- `MeteoData`: Meteorological data including atmospheric temperature, pressure, CO2 concentration
- `HumidityAtm`: Atmospheric humidity data structure with vapor pressure
- `ParVegRoof`: Roof vegetation parameters including LAI, photosynthesis properties
- `SoilPotW_ittm`: Soil water potential data for roof vegetation
- `CiCO2Leaf_ittm`: Previous iteration leaf CO2 concentration data
- `PropOpticalRoof`: Roof optical properties for radiation absorption
- `ra`: Aerodynamic resistance between roof and atmosphere [s/m]
- `rb`: Boundary layer resistance for roof vegetation [s/m]

# Returns
- `rs_sun`: Stomatal resistance for sunlit roof vegetation leaves [s/m]
- `rs_shd`: Stomatal resistance for shaded roof vegetation leaves [s/m]
- `Ci_sun`: Internal CO2 concentration for sunlit roof vegetation leaves [ppm]
- `Ci_shd`: Internal CO2 concentration for shaded roof vegetation leaves [ppm]
"""
function precalculate_stomatal_resistance_roof(
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    HumidityAtm::ModelComponents.ModelVariables.Humidity{FT},
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SoilPotW_ittm::ModelComponents.ModelVariables.SoilPotW{FT},
    CiCO2Leaf_ittm::ModelComponents.ModelVariables.CiCO2Leaf{FT},
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    ra::FT,
    rb::FT,
) where {FT<:AbstractFloat}

    # Extract input variables
    Troof_veg_tm1 = TempVec_ittm.TRoofVeg
    Tatm = MeteoData.Tatm
    Pre = MeteoData.Pre
    ea = MeteoData.ea
    Catm_O2 = MeteoData.Catm_O2
    Catm_CO2 = MeteoData.Catm_CO2
    esat_Tatm = MeteoData.AtmVapourPreSat

    # Vegetation parameters
    LAI_roof = ParVegRoof.LAI
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

    # Previous timestep values
    Psi_ltm1 = SoilPotW_ittm.SoilPotWRoofVeg_L
    Ci_sun_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegSun
    Ci_shd_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegShd

    # Parameters for stomata resistance
    Citm1_sun = Ci_sun_tm1
    Citm1_shd = Ci_shd_tm1
    Ds_atm = esat_Tatm - ea
    Oa = Catm_O2

    # Radiation calculations
    SWRabs_dir = (1 - PropOpticalRoof.aveg) * MeteoData.SW_dir
    SWRabs_diff = (1 - PropOpticalRoof.aveg) * MeteoData.SW_diff

    # Partitioning of radiation
    Fsun = (1 - exp(-Kopt_roof * LAI_roof)) / (Kopt_roof * LAI_roof)
    Fsun = clamp(Fsun, 0.0, 1.0)
    Fshd = 1 - Fsun

    PAR_sun = SWRabs_dir + Fsun * SWRabs_diff
    PAR_shd = Fshd * SWRabs_diff

    # Calculate stomatal resistance using canopy_resistance_an_evolution
    Opt_CR = 1.0  # Numerical tolerance for internal CO2 computation
    rs_sun, rs_shd, Ci_sun, Ci_shd, _, _, _, _ = canopy_resistance_an_evolution(
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
        Troof_veg_tm1 - FT(273.15),
        Tatm - FT(273.15),
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
        Opt_CR,
    )

    return rs_sun, rs_shd, Ci_sun, Ci_shd
end

function precalculate_stomatal_resistance_roof(
    TempVec_ittm::NamedTuple,
    MeteoData::NamedTuple,
    HumidityAtm::NamedTuple,
    ParVegRoof::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    SoilPotW_ittm::NamedTuple,
    CiCO2Leaf_ittm::NamedTuple,
    PropOpticalRoof::ModelComponents.Parameters.VegetatedOpticalProperties{FT},
    ra::FT,
    rb::FT,
) where {FT<:AbstractFloat}

    # Extract input variables
    Troof_veg_tm1 = TempVec_ittm.TRoofVeg
    Tatm = MeteoData.Tatm
    Pre = MeteoData.Pre
    ea = MeteoData.ea
    Catm_O2 = MeteoData.Catm_O2
    Catm_CO2 = MeteoData.Catm_CO2
    esat_Tatm = HumidityAtm.AtmVapourPreSat

    # Vegetation parameters
    LAI_roof = ParVegRoof.LAI
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

    # Previous timestep values
    Psi_ltm1 = SoilPotW_ittm.SoilPotWRoofVeg_L
    Ci_sun_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegSun
    Ci_shd_tm1 = CiCO2Leaf_ittm.CiCO2LeafRoofVegShd

    # Parameters for stomata resistance
    Citm1_sun = Ci_sun_tm1
    Citm1_shd = Ci_shd_tm1
    Ds_atm = esat_Tatm - ea
    Oa = Catm_O2

    # Radiation calculations
    SWRabs_dir = (1 - PropOpticalRoof.aveg) * MeteoData.SW_dir
    SWRabs_diff = (1 - PropOpticalRoof.aveg) * MeteoData.SW_diff

    # Partitioning of radiation
    Fsun = (1.0 - exp(-Kopt_roof * LAI_roof)) / (Kopt_roof * LAI_roof)
    Fsun = clamp(Fsun, 0.0, 1.0)
    Fshd = 1 - Fsun

    PAR_sun = SWRabs_dir + Fsun * SWRabs_diff
    PAR_shd = Fshd * SWRabs_diff

    # Calculate stomatal resistance using canopy_resistance_an_evolution
    Opt_CR = 1.0  # Numerical tolerance for internal CO2 computation
    rs_sun, rs_shd, Ci_sun, Ci_shd, _, _, _, _ = canopy_resistance_an_evolution(
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
        Troof_veg_tm1 - 273.15,
        Tatm - 273.15,
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
        Opt_CR,
    )

    return rs_sun, rs_shd, Ci_sun, Ci_shd
end
