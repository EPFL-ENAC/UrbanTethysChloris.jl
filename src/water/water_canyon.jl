"""
    water_canyon(
        MeteoData::NamedTuple,
        Int_ittm::NamedTuple,
        Owater_ittm::NamedTuple,
        Runon_ittm::NamedTuple,
        Qinlat_ittm::NamedTuple,
        Etree_In::FT,
        Egveg_In::FT,
        Egimp_Pond::FT,
        Egbare_Pond::FT,
        Egveg_Pond::FT,
        Egbare_soil::FT,
        Egveg_soil::FT,
        TEgveg::FT,
        TEtree::FT,
        ParSoilGround::NamedTuple,
        ParInterceptionTree::NamedTuple,
        ParCalculation::NamedTuple,
        ParVegGround::NamedTuple,
        ParVegTree::NamedTuple,
        FractionsGround::NamedTuple,
        Gemeotry_m::NamedTuple,
        Anthropogenic::NamedTuple
    ) where {FT<:AbstractFloat}

Calculate water balance for canyon surfaces including trees, vegetation, and impervious areas.

# Arguments
- `MeteoData`: Meteorological data
- `Int_ittm`: Previous timestep interception values
- `Owater_ittm`: Previous timestep soil moisture values
- `Runon_ittm`: Previous timestep runon values
- `Qinlat_ittm`: Previous timestep lateral flow values
- `Etree_In`: Tree evaporation [kg/m²s]
- `Egveg_In`: Ground vegetation evaporation [kg/m²s]
- `Egimp_Pond`: Impervious ground evaporation [kg/m²s]
- `Egbare_Pond`: Bare ground evaporation [kg/m²s]
- `Egveg_Pond`: Vegetated ground evaporation [kg/m²s]
- `Egbare_soil`: Bare soil evaporation [kg/m²s]
- `Egveg_soil`: Vegetated soil evaporation [kg/m²s]
- `TEgveg`: Ground vegetation transpiration [kg/m²s]
- `TEtree`: Tree transpiration [kg/m²s]
- `ParSoilGround`: Soil parameters for ground
- `ParInterceptionTree`: Tree interception parameters
- `ParCalculation`: Calculation parameters
- `ParVegGround`: Ground vegetation parameters
- `ParVegTree`: Tree vegetation parameters
- `FractionsGround`: Ground surface fractions
- `Gemeotry_m`: Geometry parameters
- `Anthropogenic`: Anthropogenic water inputs

# Returns
- `q_tree_dwn`: Tree throughfall [mm/dth]
- `In_tree`: Tree interception [mm]
- `dIn_tree_dt`: Change in tree interception [mm/dth]
- `q_gveg_dwn`: Ground vegetation throughfall [mm/dth]
- `In_gveg`: Ground vegetation interception [mm]
- `dIn_gveg_dt`: Change in ground vegetation interception [mm/dth]
- `q_gimp_runoff`: Impervious ground runoff [mm/dth]
- `In_gimp`: Impervious ground interception [mm]
- `dIn_gimp_dt`: Change in impervious ground interception [mm/dth]
- `f_inf_gimp`: Impervious ground infiltration [mm/h]
- `q_gbare_runoff`: Bare ground runoff [mm/dth]
- `In_gbare`: Bare ground interception [mm]
- `dIn_gbare_dt`: Change in bare ground interception [mm/dth]
- `f_inf_gbare`: Bare ground infiltration [mm/h]
- `q_gveg_runoff`: Vegetated ground runoff [mm/dth]
- `In_gveg_pond`: Ground vegetation ponding [mm]
- `dIn_gveg_pond_dt`: Change in ground vegetation ponding [mm/dth]
- `f_inf_gveg`: Vegetated ground infiltration [mm/h]
- `V_gimp`: Impervious ground water volume [mm]
- `O_gimp`: Impervious ground water content [-]
- `OS_gimp`: Impervious ground surface water content [-]
- `Lk_gimp`: Impervious ground leakage [mm/h]
- `Psi_s_H_gimp`: Impervious ground high soil water potential [MPa]
- `Psi_s_L_gimp`: Impervious ground low soil water potential [MPa]
- `Exwat_H_gimp`: Impervious ground high extractable water [mm]
- `Exwat_L_gimp`: Impervious ground low extractable water [mm]
- `Rd_gimp`: Impervious ground surface runoff [mm]
- `TEgveg_imp`: Ground vegetation transpiration from impervious [kg/m²s]
- `TEtree_imp`: Tree transpiration from impervious [kg/m²s]
- `Egimp_soil`: Impervious soil evaporation [kg/m²s]
- `dV_dt_gimp`: Change in impervious ground water volume [mm/dth]
- `Psi_soil_gimp`: Impervious ground soil water potential [MPa]
- `Kf_gimp`: Impervious ground hydraulic conductivity [mm/h]
- `V_gbare`: Bare ground water volume [mm]
- `O_gbare`: Bare ground water content [-]
- `OS_gbare`: Bare ground surface water content [-]
- `Lk_gbare`: Bare ground leakage [mm]
- `Psi_s_H_gbare`: Bare ground high soil water potential [MPa]
- `Psi_s_L_gbare`: Bare ground low soil water potential [MPa]
- `Exwat_H_gbare`: Bare ground high extractable water [mm]
- `Exwat_L_gbare`: Bare ground low extractable water [mm]
- `Rd_gbare`: Bare ground surface runoff [mm]
- `TEgveg_bare`: Ground vegetation transpiration from bare [kg/m²s]
- `TEtree_bare`: Tree transpiration from bare [kg/m²s]
- `Egbare_Soil`: Bare soil evaporation [kg/m²s]
- `dV_dt_gbare`: Change in bare ground water volume [mm/dth]
- `Psi_soil_gbare`: Bare ground soil water potential [MPa]
- `Kf_gbare`: Bare ground hydraulic conductivity [mm/h]
- `V_gveg`: Vegetated ground water volume [mm]
- `O_gveg`: Vegetated ground water content [-]
- `OS_gveg`: Vegetated ground surface water content [-]
- `Lk_gveg`: Vegetated ground leakage [mm]
- `Psi_s_H_gveg`: Vegetated ground high soil water potential [MPa]
- `Psi_s_L_gveg`: Vegetated ground low soil water potential [MPa]
- `Exwat_H_gveg`: Vegetated ground high extractable water [mm]
- `Exwat_L_gveg`: Vegetated ground low extractable water [mm]
- `Rd_gveg`: Vegetated ground surface runoff [mm]
- `TEgveg_veg`: Ground vegetation transpiration from vegetated [kg/m²s]
- `TEtree_veg`: Tree transpiration from vegetated [kg/m²s]
- `Egveg_Soil`: Vegetated soil evaporation [kg/m²s]
- `dV_dt_gveg`: Change in vegetated ground water volume [mm/dth]
- `Psi_soil_gveg`: Vegetated ground soil water potential [MPa]
- `Kf_gveg`: Vegetated ground hydraulic conductivity [mm/h]
- `Qin_imp`: Lateral flow to impervious ground [mm/dth]
- `Qin_bare`: Lateral flow to bare ground [mm/dth]
- `Qin_veg`: Lateral flow to vegetated ground [mm/dth]
- `Qin_bare2imp`: Lateral flow from bare to impervious [mm/dth]
- `Qin_bare2veg`: Lateral flow from bare to vegetated [mm/dth]
- `Qin_imp2bare`: Lateral flow from impervious to bare [mm/dth]
- `Qin_imp2veg`: Lateral flow from impervious to vegetated [mm/dth]
- `Qin_veg2imp`: Lateral flow from vegetated to impervious [mm/dth]
- `Qin_veg2bare`: Lateral flow from vegetated to bare [mm/dth]
- `V`: Total water volume [mm]
- `O`: Total water content [-]
- `OS`: Surface water content [-]
- `Lk`: Total leakage [mm/h]
- `Rd`: Surface runoff [mm]
- `dV_dt`: Change in total water volume [mm/dth]
- `Psi_s_L`: Low soil water potential [MPa]
- `Exwat_L`: Low extractable water [mm]
- `TEgveg_tot`: Total ground vegetation transpiration [kg/m²s]
- `Psi_s_H_tot`: Total high soil water potential [MPa]
- `Exwat_H`: Total extractable water [mm]
- `TEtree_tot`: Total tree transpiration [kg/m²s]
- `EB_TEtree`: Tree transpiration energy balance [W/m²]
- `EB_TEgveg`: Ground vegetation transpiration energy balance [W/m²]
- `WBIndv`: Individual water balance components
- `WBTot`: Total water balance components
- `Runoff`: Total runoff [mm/dth]
- `Runon`: Total runon [mm/dth]
- `Etot`: Total evapotranspiration [mm/dth]
- `DeepGLk`: Deep ground leakage [mm/dth]
- `StorageTot`: Total water storage change [mm/dth]
"""
function water_canyon(
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    Int_ittm::ModelComponents.ModelVariables.Interception{FT},
    Owater_ittm::ModelComponents.ModelVariables.Owater{FT,MR,MG},
    Runon_ittm::ModelComponents.ModelVariables.Runon{FT},
    Qinlat_ittm::ModelComponents.ModelVariables.Qinlat{FT},
    Etree_In::FT,
    Egveg_In::FT,
    Egimp_Pond::FT,
    Egbare_Pond::FT,
    Egveg_Pond::FT,
    Egbare_soil::FT,
    Egveg_soil::FT,
    TEgveg::FT,
    TEtree::FT,
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParInterceptionTree::NamedTuple,
    ParCalculation::NamedTuple,
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    Anthropogenic::ModelComponents.ForcingInputs.AnthropogenicInputs{FT,0},
) where {FT<:AbstractFloat,MR,MG}

    # Extract parameters from dictionaries
    Rain = MeteoData.Rain
    In_gimp_tm1 = Int_ittm.IntGroundImp
    In_gbare_tm1 = Int_ittm.IntGroundBare
    In_gveg_tm1 = Int_ittm.IntGroundVegPlant
    In_gvegpond_tm1 = Int_ittm.IntGroundVegGround
    In_tree_tm1 = Int_ittm.IntTree
    Otm1_imp = Owater_ittm.OwGroundSoilImp
    Otm1_bare = Owater_ittm.OwGroundSoilBare
    Otm1_veg = Owater_ittm.OwGroundSoilVeg
    Qin_imp_tm1 = Qinlat_ittm.Qin_imp
    Qin_bare_tm1 = Qinlat_ittm.Qin_bare
    Qin_veg_tm1 = Qinlat_ittm.Qin_veg
    Runon_tm1 = Runon_ittm.RunonGroundTot

    # Geometry and ground parameters
    Wcan = Gemeotry_m.Width_canyon
    LAI_g = ParVegGround.LAI
    SAI_g = ParVegGround.SAI
    LAI_t = ParVegTree.LAI
    SAI_t = ParVegTree.SAI
    r_tree = Gemeotry_m.radius_tree

    # Soil parameters
    Pcla = ParSoilGround.Pcla
    Psan = ParSoilGround.Psan
    Porg = ParSoilGround.Porg
    Kfc = ParSoilGround.Kfc
    Phy = ParSoilGround.Phy
    SPAR = ParSoilGround.SPAR
    Kbot = ParSoilGround.Kbot
    CASE_ROOT_H = ParVegTree.CASE_ROOT
    CASE_ROOT_L = ParVegGround.CASE_ROOT
    ZR95_H = ParVegTree.ZR95
    ZR95_L = ParVegGround.ZR95
    ZR50_H = ParVegTree.ZR50
    ZR50_L = ParVegGround.ZR50
    ZRmax_H = ParVegTree.ZRmax
    ZRmax_L = ParVegGround.ZRmax
    Zs = ParSoilGround.Zs
    In_max_gimp = ParSoilGround.In_max_imp
    In_max_gbare = ParSoilGround.In_max_bare
    In_max_gvegpond = ParSoilGround.In_max_underveg
    Sp_In_tree = ParInterceptionTree.Sp_In
    Sp_In_g = ParSoilGround.Sp_In
    Kimp = ParSoilGround.Kimp

    # Root parameters
    Rrootl_H = ParVegTree.Rrootl
    Rrootl_L = ParVegGround.Rrootl
    PsiL50_H = ParVegTree.PsiL50
    PsiL50_L = ParVegGround.PsiL50
    PsiX50_H = ParVegTree.PsiX50
    PsiX50_L = ParVegGround.PsiX50
    SPARTREE = ParVegTree.SPARTREE

    # Ground fractions
    Per_runoff = FractionsGround.Per_runoff
    fimp = FractionsGround.fimp
    fbare = FractionsGround.fbare
    fveg = FractionsGround.fveg

    # Time parameters
    dth = ParCalculation.dth
    row = ParCalculation.row

    # Surface presence indicators
    Cimp = FT(fimp > 0)
    Cbare = FT(fbare > 0)
    Cveg = FT(fveg > 0)
    Ctree = FT(Gemeotry_m.trees == 1)

    # Trees: water intercepted
    q_tree_dwn, In_tree, dIn_tree_dt, WB_In_tree = water_vegetation(
        Rain, Etree_In, In_tree_tm1, Sp_In_tree, LAI_t, SAI_t, row, dth
    )

    q_tree_dwn = q_tree_dwn * Ctree     # [mm/dth]
    In_tree = In_tree * Ctree           # [mm]
    dIn_tree_dt = dIn_tree_dt * Ctree   # [mm/dth]
    WB_In_tree = WB_In_tree * Ctree     # [mm/dth]

    # Water received by any ground fraction including rain and dripping from trees
    Rain_ground = 4 * r_tree * Ctree * q_tree_dwn + (1 - 4 * r_tree * Ctree) * Rain  # [mm/dth]

    # Ground vegetation: interception
    q_gveg_dwn, In_gveg, dIn_gveg_dt, WB_In_gveg = water_vegetation(
        Rain_ground, Egveg_In, In_gveg_tm1, Sp_In_g, LAI_g, SAI_g, row, dth
    )

    q_gveg_dwn = q_gveg_dwn * Cveg      # [mm/dth]
    In_gveg = In_gveg * Cveg            # [mm]
    dIn_gveg_dt = dIn_gveg_dt * Cveg    # [mm/dth]
    WB_In_gveg = WB_In_gveg * Cveg      # [mm/dth]

    # Ground water calculations
    f_inf_gimp = Kimp  # [mm/h]

    # Impervious ground
    q_gimp_runoff, In_gimp, dIn_gimp_dt, f_inf_gimp, WB_In_gimp = water_impervious(
        Rain_ground, Runon_tm1, Egimp_Pond, In_gimp_tm1, dth, row, In_max_gimp, f_inf_gimp
    )

    q_gimp_runoff = q_gimp_runoff * Cimp    # [mm/dth]
    In_gimp = In_gimp * Cimp                # [mm/dth]
    dIn_gimp_dt = dIn_gimp_dt * Cimp        # [mm/dth]
    f_inf_gimp = f_inf_gimp * Cimp          # [mm/h]
    WB_In_gimp = WB_In_gimp * Cimp          # [mm/dth]

    # Bare ground
    q_gbare_runoff, In_gbare, dIn_gbare_dt, f_inf_gbare, WB_In_gbare = water_ground(
        Rain_ground + Anthropogenic.Waterf_canyonBare,
        Runon_tm1,
        Egbare_Pond,
        Otm1_bare,
        In_gbare_tm1,
        In_max_gbare,
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
        ZR95_L .* NaN,
        ZR50_H,
        ZR50_L .* NaN,
        ZRmax_H,
        ZRmax_L .* NaN,
        Zs,
        dth,
        row,
    )

    q_gbare_runoff = q_gbare_runoff * Cbare    # [mm/dth]
    In_gbare = In_gbare * Cbare                # [mm/dth]
    dIn_gbare_dt = dIn_gbare_dt * Cbare        # [mm/dth]
    f_inf_gbare = f_inf_gbare * Cbare          # [mm/h]
    WB_In_gbare = WB_In_gbare * Cbare          # [mm/dth]

    # Vegetated ground
    q_gveg_runoff, In_gveg_pond, dIn_gveg_pond_dt, f_inf_gveg, WB_Pond_gveg = water_ground(
        q_gveg_dwn + Anthropogenic.Waterf_canyonVeg,
        Runon_tm1,
        Egveg_Pond,
        Otm1_veg,
        In_gvegpond_tm1,
        In_max_gvegpond,
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
        dth,
        row,
    )

    q_gveg_runoff = q_gveg_runoff * Cveg          # [mm/dth]
    In_gveg_pond = In_gveg_pond * Cveg            # [mm/dth]
    dIn_gveg_pond_dt = dIn_gveg_pond_dt * Cveg    # [mm/dth]
    f_inf_gveg = f_inf_gveg * Cveg                # [mm/h]
    WB_Pond_gveg = WB_Pond_gveg * Cveg            # [mm/dth]

    # Soil water redistribution
    if SPARTREE == 1
        TEtree0_imp = TEtree * Ctree * (4 * r_tree) * Cimp
        TEtree0_bare = TEtree * Ctree * (4 * r_tree) * Cbare
        TEtree0_veg = TEtree * Ctree * (4 * r_tree) * Cveg
    else
        if (4 * r_tree) <= (fveg + fbare)
            TEtree0_imp = zero(FT) * Ctree * Cimp
            TEtree0_bare = TEtree * (4 * r_tree) / (fveg + fbare) * Ctree * Cbare
            TEtree0_veg = TEtree * (4 * r_tree) / (fveg + fbare) * Ctree * Cveg
        else
            TEtree0_imp = ((4 * r_tree) - (fveg + fbare)) * TEtree / fimp * Ctree * Cimp
            TEtree0_bare = TEtree * Ctree * Cbare
            TEtree0_veg = TEtree * Ctree * Cveg
        end
    end

    TEgveg = TEgveg * Cveg
    Egbare_soil = Egbare_soil * Cbare
    Egveg_soil = Egveg_soil * Cveg

    # Impervious soil column (from layer 3)
    V_gimp1, O_gimp1, _, Lk_gimp1, _, _, Exwat_H_gimp1, Exwat_L_gimp1, Rd_gimp1, TEgveg_imp1, TEtree_imp1, Egimp_soil1, _, WB_Soil_gimp1, Psi_Soil_gimp1, Kf_gimp1 = water_soil(
        Otm1_imp[3:end],
        f_inf_gimp,
        TEtree0_imp,
        zero(FT),
        zero(FT),
        zeros(FT, length(Qin_imp_tm1[3:end])),
        dth,
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
        [zero(FT)],
        Rrootl_H,
        [zero(FT)],
        PsiL50_H,
        [zero(FT)],
        PsiX50_H,
        Zs[3:end] .- Zs[3],
        row,
    )

    V_gimp1 = [NaN; NaN; V_gimp1]

    # Bare soil column
    V_gbare1, _, _, Lk_gbare1, _, _, _, _, Rd_gbare1, TEgveg_bare1, TEtree_bare1, Egbare_Soil1, _, WB_Soil_gbare1, _, _ = water_soil(
        Otm1_bare,
        f_inf_gbare,
        TEtree0_bare,
        zero(FT),
        Egbare_soil,
        zeros(FT, length(Qin_bare_tm1)),
        dth,
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
        [zero(FT)],
        Rrootl_H,
        [zero(FT)],
        PsiL50_H,
        [zero(FT)],
        PsiX50_H,
        Zs,
        row,
    )

    # Vegetated soil column
    V_gveg1, _, _, Lk_gveg1, _, _, _, _, Rd_gveg1, TEgveg_veg1, TEtree_veg1, Egveg_Soil1, _, WB_Soil_gveg1, _, _ = water_soil(
        Otm1_veg,
        f_inf_gveg,
        TEtree0_veg,
        TEgveg,
        Egveg_soil,
        zeros(FT, length(Qin_veg_tm1)),
        dth,
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
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
        Zs,
        row,
    )

    # Lateral water redistribution parameters
    T_SPAN = (zero(FT), dth)
    OPT_SM = (abstol=FT(0.05), maxstep=dth)
    dz = diff(Zs)

    # Get soil parameters
    _, _, _, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, _, _, _, _, _, _, _, _, _, _ = Soil.soil_parameters_total(
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

    # Initialize redistributed volumes
    V_gimp2 = zeros(FT, length(dz))
    V_gbare2 = zeros(FT, length(dz))
    V_gveg2 = zeros(FT, length(dz))

    # First two layers: exchange between bare and vegetated only
    for i in 1:2
        Vlat1 = [V_gbare1[i], V_gveg1[i]]

        sol = solve(
            ODEProblem(
                (V, p, t) -> Soil.soil_moistures_rich_comp_lat2(
                    V,
                    dz[i],
                    SPAR,
                    Ks_Zs[i],
                    Osat[i],
                    Ohy[i],
                    L[i],
                    Pe[i],
                    O33[i],
                    alpVG[i],
                    nVG[i],
                    Cbare,
                    Cveg,
                    fbare,
                    fveg,
                    Wcan,
                ),
                Vlat1,
                T_SPAN,
            ),
            Rosenbrock23(; autodiff=AutoFiniteDiff());
            abstol=OPT_SM.abstol,
            dtmax=OPT_SM.maxstep,
        )

        V_gimp2[i] = NaN
        V_gbare2[i] = sol.u[end][1]
        V_gveg2[i] = sol.u[end][2]
    end

    # Remaining layers: exchange between all three columns
    for i in 3:length(dz)
        Vlat1 = [V_gimp1[i], V_gbare1[i], V_gveg1[i]]

        sol = solve(
            ODEProblem(
                (V, p, t) -> Soil.soil_moistures_rich_comp_lat3(
                    V,
                    dz[i],
                    SPAR,
                    Ks_Zs[i],
                    Osat[i],
                    Ohy[i],
                    L[i],
                    Pe[i],
                    O33[i],
                    alpVG[i],
                    nVG[i],
                    Cimp,
                    Cbare,
                    Cveg,
                    fimp,
                    fbare,
                    fveg,
                    Wcan,
                ),
                Vlat1,
                T_SPAN,
            ),
            Rosenbrock23(; autodiff=AutoFiniteDiff());
            abstol=OPT_SM.abstol,
            dtmax=OPT_SM.maxstep,
        )

        V_gimp2[i] = sol.u[end][1]
        V_gbare2[i] = sol.u[end][2]
        V_gveg2[i] = sol.u[end][3]
    end

    # Back compute lateral fluxes
    Qin_imp = V_gimp2 .- V_gimp1            # [mm/dth]
    Qin_bare = V_gbare2 .- V_gbare1         # [mm/dth]
    Qin_veg = V_gveg2 .- V_gveg1           # [mm/dth]

    Qin_bare2imp = zeros(FT, length(dz))
    Qin_bare2veg = zeros(FT, length(dz))
    Qin_imp2bare = zeros(FT, length(dz))
    Qin_imp2veg = zeros(FT, length(dz))
    Qin_veg2imp = zeros(FT, length(dz))
    Qin_veg2bare = zeros(FT, length(dz))

    Qin_bare2imp[1:2] .= NaN
    Qin_imp2bare[1:2] .= NaN
    Qin_imp2veg[1:2] .= NaN
    Qin_veg2imp[1:2] .= NaN

    # Impervious
    V_gimp2, O_gimp2, OS_gimp2, Psi_Soil_gimp2, Psi_s_H_gimp2, Psi_s_L_gimp2, Exwat_H_gimp2, Exwat_L_gimp2, Kf_gimp2 = Soil.soil_moisture_conductivity_update(
        V_gimp2[3:end],
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
        Zs[3:end] .- Zs[3],
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    V_gimp2 = [NaN; NaN; V_gimp2]
    O_gimp2 = [NaN; NaN; O_gimp2]
    Exwat_H_gimp2 = [NaN; NaN; vec(Exwat_H_gimp2)]
    Exwat_L_gimp2 = [NaN; NaN; vec(Exwat_L_gimp2)]
    Psi_Soil_gimp2 = [NaN; NaN; Psi_Soil_gimp2]
    Kf_gimp2 = [NaN; NaN; Kf_gimp2]

    # Bare
    V_gbare2, O_gbare2, OS_gbare2, Psi_soil_gbare2, Psi_s_H_gbare2, Psi_s_L_gbare2, Exwat_H_gbare2, Exwat_L_gbare2, Kf_gbare2 = Soil.soil_moisture_conductivity_update(
        V_gbare2,
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
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    # Vegetated
    V_gveg2, O_gveg2, OS_gveg2, Psi_soil_gveg2, Psi_s_H_gveg2, Psi_s_L_gveg2, Exwat_H_gveg2, Exwat_L_gveg2, Kf_gveg2 = Soil.soil_moisture_conductivity_update(
        V_gveg2,
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
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    # Change in water volume
    Vtm1_imp = (Otm1_imp - Ohy) .* dz
    Vtm1_imp[1:2] .= NaN
    Vtm1_bare = (Otm1_bare - Ohy) .* dz
    Vtm1_veg = (Otm1_veg - Ohy) .* dz

    function nansum(x)
        return sum(filter(!isnan, x))
    end

    function nansum(arr, dims)
        return sum(x -> !isnan(x) * x, arr; dims=dims)
    end

    dV_dt_gimpTot = nansum(V_gimp2) - nansum(Vtm1_imp)
    dV_dt_gbareTot = nansum(V_gbare2) - nansum(Vtm1_bare)
    dV_dt_gvegTot = nansum(V_gveg2) - nansum(Vtm1_veg)

    # Surface water balance: Tree, impervious, bare, vegetated ground
    WBsurf_tree =
        Ctree * Rain - Etree_In * dth * 3600 * 1000 / row - q_tree_dwn - dIn_tree_dt # Is not used [mm/dth]

    WBsurf_imp =
        Cimp * Rain_ground + Cimp * Runon_tm1 - Egimp_Pond * dth * 3600 * 1000 / row -
        q_gimp_runoff - f_inf_gimp * dth - dIn_gimp_dt # [mm/dth]

    WBsurf_bare =
        Cbare * Rain_ground + Cbare * Runon_tm1 + Anthropogenic.Waterf_canyonBare -
        Egbare_Pond * dth * 3600 * 1000 / row - f_inf_gbare * dth - q_gbare_runoff -
        dIn_gbare_dt # [mm/dth]

    WBsurf_veg =
        Cveg * Rain_ground + Cveg * Runon_tm1 + Anthropogenic.Waterf_canyonVeg -
        (Egveg_In + Egveg_Pond) * dth * 3600 * 1000 / row - f_inf_gveg * dth -
        q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt # [mm/dth]

    # Soil water balance: impervious, bare, vegetated ground
    WBsoil_imp =
        f_inf_gimp * dth + nansum(Qin_imp) .-
        (TEgveg_imp1 .+ TEtree_imp1 .+ Egimp_soil1) * dth * 3600 * 1000 / row -
        Lk_gimp1 * dth - Rd_gimp1 - dV_dt_gimpTot # [mm/dth]

    WBsoil_bare =
        f_inf_gbare * dth + nansum(Qin_bare) .-
        (TEgveg_bare1 .+ TEtree_bare1 + Egbare_Soil1) * dth * 3600 * 1000 / row -
        Lk_gbare1 * dth - Rd_gbare1 - dV_dt_gbareTot # [mm/dth]

    WBsoil_veg =
        f_inf_gveg * dth + nansum(Qin_veg) .-
        (TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1) * dth * 3600 * 1000 / row -
        Lk_gveg1 * dth - Rd_gveg1 - dV_dt_gvegTot # [mm/dth]

    # Ground water balance: impervious, bare, vegetated ground
    WBimp_tot =
        Cimp * Rain_ground + Cimp * Runon_tm1 + nansum(Qin_imp) -
        (Egimp_Pond + TEgveg_imp1 + TEtree_imp1 + Egimp_soil1) * dth * 3600 * 1000 / row -
        q_gimp_runoff - dIn_gimp_dt - Lk_gimp1 * dth - Rd_gimp1 - dV_dt_gimpTot # [mm/dth]

    WBbare_tot =
        Cbare * Rain_ground +
        Cbare * Runon_tm1 +
        nansum(Qin_bare) +
        Anthropogenic.Waterf_canyonBare -
        (Egbare_Pond + TEgveg_bare1 + TEtree_bare1 + Egbare_Soil1) * dth * 3600 * 1000 /
        row - q_gbare_runoff - dIn_gbare_dt - Lk_gbare1 * dth - Rd_gbare1 - dV_dt_gbareTot # [mm/dth]

    WBveg_tot =
        Cveg * Rain_ground +
        Cveg * Runon_tm1 +
        nansum(Qin_veg) +
        Anthropogenic.Waterf_canyonVeg -
        (Egveg_In + Egveg_Pond + TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1) *
        dth *
        3600 *
        1000 / row - q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt - Lk_gveg1 * dth -
        Rd_gveg1 - dV_dt_gvegTot # [mm/dth]

    # Calculate total runoff including saturated excess
    Runoff =
        Per_runoff * (
            fimp * (q_gimp_runoff + Rd_gimp1) +
            fbare * (q_gbare_runoff + Rd_gbare1) +
            fveg * (q_gveg_runoff + Rd_gveg1)
        )  # [mm/dth]

    # Calculate total runon
    Runon =
        (1 - Per_runoff) * (
            fimp * (q_gimp_runoff + Rd_gimp1) +
            fbare * (q_gbare_runoff + Rd_gbare1) +
            fveg * (q_gveg_runoff + Rd_gveg1)
        )  # [mm/dth]

    # Calculate total evapotranspiration and leakage
    Etot =
        (
            fimp * Egimp_Pond +
            fbare * (Egbare_Pond + Egbare_soil) +
            fveg * (Egveg_Pond + Egveg_soil + TEgveg) +
            TEtree
        ) *
        3600 *
        1000 / row  # [mm/h]

    # Deep ground leakage and change in storage
    DeepGLk = fimp * Lk_gimp1 + fbare * Lk_gbare1 + fveg * Lk_gveg1  # [mm/h]

    StorageTot =
        fimp * (dIn_gimp_dt + dV_dt_gimpTot) +
        fbare * (dIn_gbare_dt + dV_dt_gbareTot) +
        fveg * (dIn_gveg_dt + dIn_gveg_pond_dt + dV_dt_gvegTot) +
        (4 * r_tree) * dIn_tree_dt  # [mm/dth]

    WBcanyon_flux =
        Rain +
        Runon_tm1 +
        fveg * Anthropogenic.Waterf_canyonVeg +
        fbare * Anthropogenic.Waterf_canyonBare - Etot - DeepGLk - Runoff - Runon -
        StorageTot  # [mm]

    WBtree_level =
        (4 * r_tree) *
        (Ctree * Rain - Etree_In * dth * 3600 * 1000 / row - q_tree_dwn - dIn_tree_dt)  # [mm/dth]

    WBground_level =
        fimp * (
            Rain_ground + Runon_tm1 - dIn_gimp_dt - Egimp_Pond * dth * 3600 * 1000 / row -
            q_gimp_runoff - Rd_gimp1 - f_inf_gimp * dth
        ) +
        fbare * (
            Rain_ground + Runon_tm1 + Anthropogenic.Waterf_canyonBare - dIn_gbare_dt -
            Egbare_Pond * dth * 3600 * 1000 / row - q_gbare_runoff - Rd_gbare1 -
            f_inf_gbare * dth
        ) +
        fveg * (
            Rain_ground + Runon_tm1 + Anthropogenic.Waterf_canyonVeg - dIn_gveg_dt -
            dIn_gveg_pond_dt - (Egveg_In + Egveg_Pond) * dth * 3600 * 1000 / row -
            q_gveg_runoff - Rd_gveg1 - f_inf_gveg * dth
        )  # [mm/dth]

    WBsoil_level =
        fimp * (
            f_inf_gimp * dth +
            3600 * 1000 / row * dth * (-TEgveg_imp1 - TEtree_imp1 - Egimp_soil1) -
            Lk_gimp1 * dth - dV_dt_gimpTot
        ) +
        fbare * (
            f_inf_gbare * dth +
            3600 * 1000 / row * dth * (-TEgveg_bare1 - TEtree_bare1 - Egbare_Soil1) -
            Lk_gbare1 * dth - dV_dt_gbareTot
        ) +
        fveg * (
            f_inf_gveg * dth +
            3600 * 1000 / row * dth * (-TEgveg_veg1 - TEtree_veg1 - Egveg_Soil1) -
            Lk_gveg1 * dth - dV_dt_gvegTot
        )  # [mm/dth]

    WBcanyon_level =
        Rain +
        Runon_tm1 +
        (4 * r_tree) * (-Etree_In * dth * 3600 * 1000 / row - dIn_tree_dt) +
        fimp * (
            -dIn_gimp_dt - Egimp_Pond * dth * 3600 * 1000 / row - q_gimp_runoff - Rd_gimp1 -
            3600 * 1000 * dth / row * (TEgveg_imp1 + TEtree_imp1 + Egimp_soil1) -
            Lk_gimp1 * dth - dV_dt_gimpTot
        ) +
        fbare * (
            Anthropogenic.Waterf_canyonBare - dIn_gbare_dt -
            Egbare_Pond * dth * 3600 * 1000 / row - q_gbare_runoff - Rd_gbare1 -
            3600 * 1000 / row * dth * (TEgveg_bare1 + TEtree_bare1 + Egbare_Soil1) -
            Lk_gbare1 * dth - dV_dt_gbareTot
        ) +
        fveg * (
            Anthropogenic.Waterf_canyonVeg - dIn_gveg_dt - dIn_gveg_pond_dt -
            (Egveg_In + Egveg_Pond) * dth * 3600 * 1000 / row - q_gveg_runoff - Rd_gveg1 -
            3600 * 1000 / row * dth * (TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1) -
            Lk_gveg1 * dth - dV_dt_gvegTot
        )  # [mm/dth]

    # Assign final variables
    V_gimp = V_gimp2
    O_gimp = O_gimp2
    V_gimp[isnan.(V_gimp)] .= 0
    O_gimp[isnan.(O_gimp)] .= 0
    OS_gimp = OS_gimp2
    Lk_gimp = Lk_gimp1
    Psi_s_H_gimp = Psi_s_H_gimp2
    Psi_s_L_gimp = Psi_s_L_gimp2
    Exwat_H_gimp = vec(Exwat_H_gimp2)
    Exwat_L_gimp = vec(Exwat_L_gimp2)
    Exwat_H_gimp[isnan.(Exwat_H_gimp)] .= 0
    Exwat_L_gimp[isnan.(Exwat_L_gimp)] .= 0
    Rd_gimp = Rd_gimp1
    TEgveg_imp = TEgveg_imp1
    Egimp_soil = Egimp_soil1
    dV_dt_gimp = dV_dt_gimpTot
    Psi_soil_gimp = Psi_Soil_gimp2
    Kf_gimp = Kf_gimp2

    # Update bare soil variables
    V_gbare = V_gbare2
    O_gbare = O_gbare2
    OS_gbare = OS_gbare2
    Lk_gbare = Lk_gbare1
    Psi_s_H_gbare = Psi_s_H_gbare2
    Psi_s_L_gbare = Psi_s_L_gbare2
    Exwat_H_gbare = Exwat_H_gbare2
    Exwat_L_gbare = Exwat_L_gbare2
    Rd_gbare = Rd_gbare1
    TEgveg_bare = TEgveg_bare1
    Egbare_Soil = Egbare_Soil1
    dV_dt_gbare = dV_dt_gbareTot
    Psi_soil_gbare = Psi_soil_gbare2
    Kf_gbare = Kf_gbare2

    # Update vegetated soil variables
    V_gveg = V_gveg2
    O_gveg = O_gveg2
    OS_gveg = OS_gveg2
    Lk_gveg = Lk_gveg1
    Psi_s_H_gveg = Psi_s_H_gveg2
    Psi_s_L_gveg = Psi_s_L_gveg2
    Exwat_H_gveg = Exwat_H_gveg2
    Exwat_L_gveg = Exwat_L_gveg2
    Rd_gveg = Rd_gveg1
    TEgveg_veg = TEgveg_veg1
    Egveg_Soil = Egveg_Soil1
    dV_dt_gveg = dV_dt_gvegTot
    Psi_soil_gveg = Psi_soil_gveg2
    Kf_gveg = Kf_gveg2

    # Update water balance variables
    WB_Soil_gimp = WB_Soil_gimp1
    WB_Soil_gbare = WB_Soil_gbare1
    WB_Soil_gveg = WB_Soil_gveg1

    # Update total values and rescale
    V = nansum(hcat(fimp * V_gimp, fbare * V_gbare, fveg * V_gveg), 2)
    O = nansum(hcat(fimp * O_gimp, fbare * O_gbare, fveg * O_gveg), 2)
    OS = nansum(hcat(fimp * OS_gimp, fbare * OS_gbare, fveg * OS_gveg), 2)
    Lk = nansum(hcat(fimp * Lk_gimp, fbare * Lk_gbare, fveg * Lk_gveg), 2)
    Rd = nansum(hcat(fimp * Rd_gimp, fbare * Rd_gbare, fveg * Rd_gveg), 2)
    dV_dt = nansum(hcat(fimp * dV_dt_gimp, fbare * dV_dt_gbare, fveg * dV_dt_gveg), 2)

    # Rescale first two layers
    V[1:2] ./= (fbare + fveg)
    O[1:2] ./= (fbare + fveg)
    V[isnan.(V)] .= 0
    O[isnan.(O)] .= 0

    # Assign vegetation parameters
    Psi_s_L = Psi_s_L_gveg
    Exwat_L = Exwat_L_gveg
    TEgveg_tot = TEgveg_veg

    Psi_s_H_tot = Psi_s_H_gveg
    Exwat_H = @. fimp * Exwat_H_gimp + fbare * Exwat_H_gbare + fveg * Exwat_H_gveg

    # Handle tree transpiration based on SPARTREE case
    if SPARTREE == 1  # Tree roots can access all water in the soil
        TEtree_imp = TEtree_imp1 / (4 * r_tree) * Cimp
        TEtree_bare = TEtree_bare1 / (4 * r_tree) * Cbare
        TEtree_veg = TEtree_veg1 / (4 * r_tree) * Cveg
    else  # Tree roots access based on crown size
        if (4 * r_tree) <= (fveg + fbare)
            TEtree_imp = 0 * Cimp
            TEtree_bare = TEtree_bare1 / (4 * r_tree) / (fveg + fbare) * Cbare
            TEtree_veg = TEtree_veg1 / (4 * r_tree) / (fveg + fbare) * Cveg
        else
            TEtree_imp = TEtree_imp1 / ((4 * r_tree) - (fveg + fbare)) * fimp * Cimp
            TEtree_bare = TEtree_bare1 * Cbare
            TEtree_veg = TEtree_veg1 * Cveg
        end
    end

    # Calculate total tree transpiration
    TEtree_tot = fimp * TEtree_imp + fbare * TEtree_bare + fveg * TEtree_veg

    # Compute energy balance terms
    EB_TEtree = TEtree_tot - TEtree
    EB_TEgveg = TEgveg_tot - TEgveg

    # Create water balance structs using NamedTuples
    WBIndv = (
        WB_In_tree=WB_In_tree,
        WB_In_gveg=WB_In_gveg,
        WB_In_gimp=WB_In_gimp,
        WB_In_gbare=WB_In_gbare,
        WB_Pond_gveg=WB_Pond_gveg,
        WB_Soil_gimp=WB_Soil_gimp,
        WB_Soil_gbare=WB_Soil_gbare,
        WB_Soil_gveg=WB_Soil_gveg,
    )

    WBTot = (
        WBsurf_tree=WBsurf_tree,
        WBsurf_imp=WBsurf_imp,
        WBsurf_bare=WBsurf_bare,
        WBsurf_veg=WBsurf_veg,
        WBsoil_imp=WBsoil_imp,
        WBsoil_bare=WBsoil_bare,
        WBsoil_veg=WBsoil_veg,
        WBimp_tot=WBimp_tot,
        WBbare_tot=WBbare_tot,
        WBveg_tot=WBveg_tot,
        WBcanyon_flux=WBcanyon_flux,
        WBtree_level=WBtree_level,
        WBground_level=WBground_level,
        WBsoil_level=WBsoil_level,
        WBcanyon_level=WBcanyon_level,
    )

    return (;
        q_tree_dwn,
        In_tree,
        dIn_tree_dt,
        q_gveg_dwn,
        In_gveg,
        dIn_gveg_dt,
        q_gimp_runoff,
        In_gimp,
        dIn_gimp_dt,
        f_inf_gimp,
        q_gbare_runoff,
        In_gbare,
        dIn_gbare_dt,
        f_inf_gbare,
        q_gveg_runoff,
        In_gveg_pond,
        dIn_gveg_pond_dt,
        f_inf_gveg,
        V_gimp,
        O_gimp,
        OS_gimp,
        Lk_gimp,
        Psi_s_H_gimp,
        Psi_s_L_gimp,
        Exwat_H_gimp,
        Exwat_L_gimp,
        Rd_gimp,
        TEgveg_imp,
        TEtree_imp,
        Egimp_soil,
        dV_dt_gimp,
        Psi_soil_gimp,
        Kf_gimp,
        V_gbare,
        O_gbare,
        OS_gbare,
        Lk_gbare,
        Psi_s_H_gbare,
        Psi_s_L_gbare,
        Exwat_H_gbare,
        Exwat_L_gbare,
        Rd_gbare,
        TEgveg_bare,
        TEtree_bare,
        Egbare_Soil,
        dV_dt_gbare,
        Psi_soil_gbare,
        Kf_gbare,
        V_gveg,
        O_gveg,
        OS_gveg,
        Lk_gveg,
        Psi_s_H_gveg,
        Psi_s_L_gveg,
        Exwat_H_gveg,
        Exwat_L_gveg,
        Rd_gveg,
        TEgveg_veg,
        TEtree_veg,
        Egveg_Soil,
        dV_dt_gveg,
        Psi_soil_gveg,
        Kf_gveg,
        Qin_imp,
        Qin_bare,
        Qin_veg,
        Qin_bare2imp,
        Qin_bare2veg,
        Qin_imp2bare,
        Qin_imp2veg,
        Qin_veg2imp,
        Qin_veg2bare,
        V,
        O,
        OS,
        Lk,
        Rd,
        dV_dt,
        Psi_s_L,
        Exwat_L,
        TEgveg_tot,
        Psi_s_H_tot,
        Exwat_H,
        TEtree_tot,
        EB_TEtree,
        EB_TEgveg,
        WBIndv,
        WBTot,
        Runoff,
        Runon_ittm,
        Etot,
        DeepGLk,
        StorageTot,
    )
end

function water_canyon(
    MeteoData::NamedTuple,
    Int_ittm::NamedTuple,
    Owater_ittm::NamedTuple,
    Runon_ittm::NamedTuple,
    Qinlat_ittm::NamedTuple,
    Etree_In::FT,
    Egveg_In::FT,
    Egimp_Pond::FT,
    Egbare_Pond::FT,
    Egveg_Pond::FT,
    Egbare_soil::FT,
    Egveg_soil::FT,
    TEgveg::FT,
    TEtree::FT,
    ParSoilGround::ModelComponents.Parameters.VegetatedSoilParameters{FT},
    ParInterceptionTree::NamedTuple,
    ParCalculation::NamedTuple,
    ParVegGround::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    ParVegTree::ModelComponents.Parameters.HeightDependentVegetationParameters{FT},
    FractionsGround::ModelComponents.Parameters.LocationSpecificSurfaceFractions{FT},
    Gemeotry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    Anthropogenic::NamedTuple,
) where {FT<:AbstractFloat}

    # Extract parameters from dictionaries
    Rain = MeteoData.Rain
    In_gimp_tm1 = Int_ittm.IntGroundImp
    In_gbare_tm1 = Int_ittm.IntGroundBare
    In_gveg_tm1 = Int_ittm.IntGroundVegPlant
    In_gvegpond_tm1 = Int_ittm.IntGroundVegGround
    In_tree_tm1 = Int_ittm.IntTree
    Otm1_imp = Owater_ittm.OwGroundSoilImp
    Otm1_bare = Owater_ittm.OwGroundSoilBare
    Otm1_veg = Owater_ittm.OwGroundSoilVeg
    Qin_imp_tm1 = Qinlat_ittm.Qin_imp
    Qin_bare_tm1 = Qinlat_ittm.Qin_bare
    Qin_veg_tm1 = Qinlat_ittm.Qin_veg
    Runon_tm1 = Runon_ittm.RunonGroundTot

    # Geometry and ground parameters
    Wcan = Gemeotry_m.Width_canyon
    LAI_g = ParVegGround.LAI
    SAI_g = ParVegGround.SAI
    LAI_t = ParVegTree.LAI
    SAI_t = ParVegTree.SAI
    r_tree = Gemeotry_m.radius_tree

    # Soil parameters
    Pcla = ParSoilGround.Pcla
    Psan = ParSoilGround.Psan
    Porg = ParSoilGround.Porg
    Kfc = ParSoilGround.Kfc
    Phy = ParSoilGround.Phy
    SPAR = ParSoilGround.SPAR
    Kbot = ParSoilGround.Kbot
    CASE_ROOT_H = ParVegTree.CASE_ROOT
    CASE_ROOT_L = ParVegGround.CASE_ROOT
    ZR95_H = ParVegTree.ZR95
    ZR95_L = ParVegGround.ZR95
    ZR50_H = ParVegTree.ZR50
    ZR50_L = ParVegGround.ZR50
    ZRmax_H = ParVegTree.ZRmax
    ZRmax_L = ParVegGround.ZRmax
    Zs = ParSoilGround.Zs
    In_max_gimp = ParSoilGround.In_max_imp
    In_max_gbare = ParSoilGround.In_max_bare
    In_max_gvegpond = ParSoilGround.In_max_underveg
    Sp_In_tree = ParInterceptionTree.Sp_In
    Sp_In_g = ParSoilGround.Sp_In
    Kimp = ParSoilGround.Kimp

    # Root parameters
    Rrootl_H = ParVegTree.Rrootl
    Rrootl_L = ParVegGround.Rrootl
    PsiL50_H = ParVegTree.PsiL50
    PsiL50_L = ParVegGround.PsiL50
    PsiX50_H = ParVegTree.PsiX50
    PsiX50_L = ParVegGround.PsiX50
    SPARTREE = ParVegTree.SPARTREE

    # Ground fractions
    Per_runoff = FractionsGround.Per_runoff
    fimp = FractionsGround.fimp
    fbare = FractionsGround.fbare
    fveg = FractionsGround.fveg

    # Time parameters
    dth = ParCalculation.dth
    row = ParCalculation.row

    # Surface presence indicators
    Cimp = FT(fimp > 0)
    Cbare = FT(fbare > 0)
    Cveg = FT(fveg > 0)
    Ctree = FT(Gemeotry_m.trees == 1)

    # Trees: water intercepted
    q_tree_dwn, In_tree, dIn_tree_dt, WB_In_tree = water_vegetation(
        Rain, Etree_In, In_tree_tm1, Sp_In_tree, LAI_t, SAI_t, row, dth
    )

    q_tree_dwn = q_tree_dwn * Ctree     # [mm/dth]
    In_tree = In_tree * Ctree           # [mm]
    dIn_tree_dt = dIn_tree_dt * Ctree   # [mm/dth]
    WB_In_tree = WB_In_tree * Ctree     # [mm/dth]

    # Water received by any ground fraction including rain and dripping from trees
    Rain_ground = 4 * r_tree * Ctree * q_tree_dwn + (1 - 4 * r_tree * Ctree) * Rain  # [mm/dth]

    # Ground vegetation: interception
    q_gveg_dwn, In_gveg, dIn_gveg_dt, WB_In_gveg = water_vegetation(
        Rain_ground, Egveg_In, In_gveg_tm1, Sp_In_g, LAI_g, SAI_g, row, dth
    )

    q_gveg_dwn = q_gveg_dwn * Cveg      # [mm/dth]
    In_gveg = In_gveg * Cveg            # [mm]
    dIn_gveg_dt = dIn_gveg_dt * Cveg    # [mm/dth]
    WB_In_gveg = WB_In_gveg * Cveg      # [mm/dth]

    # Ground water calculations
    f_inf_gimp = Kimp  # [mm/h]

    # Impervious ground
    q_gimp_runoff, In_gimp, dIn_gimp_dt, f_inf_gimp, WB_In_gimp = water_impervious(
        Rain_ground, Runon_tm1, Egimp_Pond, In_gimp_tm1, dth, row, In_max_gimp, f_inf_gimp
    )

    q_gimp_runoff = q_gimp_runoff * Cimp    # [mm/dth]
    In_gimp = In_gimp * Cimp                # [mm/dth]
    dIn_gimp_dt = dIn_gimp_dt * Cimp        # [mm/dth]
    f_inf_gimp = f_inf_gimp * Cimp          # [mm/h]
    WB_In_gimp = WB_In_gimp * Cimp          # [mm/dth]

    # Bare ground
    q_gbare_runoff, In_gbare, dIn_gbare_dt, f_inf_gbare, WB_In_gbare = water_ground(
        Rain_ground + Anthropogenic.Waterf_canyonBare,
        Runon_tm1,
        Egbare_Pond,
        Otm1_bare,
        In_gbare_tm1,
        In_max_gbare,
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
        ZR95_L .* NaN,
        ZR50_H,
        ZR50_L .* NaN,
        ZRmax_H,
        ZRmax_L .* NaN,
        Zs,
        dth,
        row,
    )

    q_gbare_runoff = q_gbare_runoff * Cbare    # [mm/dth]
    In_gbare = In_gbare * Cbare                # [mm/dth]
    dIn_gbare_dt = dIn_gbare_dt * Cbare        # [mm/dth]
    f_inf_gbare = f_inf_gbare * Cbare          # [mm/h]
    WB_In_gbare = WB_In_gbare * Cbare          # [mm/dth]

    # Vegetated ground
    q_gveg_runoff, In_gveg_pond, dIn_gveg_pond_dt, f_inf_gveg, WB_Pond_gveg = water_ground(
        q_gveg_dwn + Anthropogenic.Waterf_canyonVeg,
        Runon_tm1,
        Egveg_Pond,
        Otm1_veg,
        In_gvegpond_tm1,
        In_max_gvegpond,
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
        dth,
        row,
    )

    q_gveg_runoff = q_gveg_runoff * Cveg          # [mm/dth]
    In_gveg_pond = In_gveg_pond * Cveg            # [mm/dth]
    dIn_gveg_pond_dt = dIn_gveg_pond_dt * Cveg    # [mm/dth]
    f_inf_gveg = f_inf_gveg * Cveg                # [mm/h]
    WB_Pond_gveg = WB_Pond_gveg * Cveg            # [mm/dth]

    # Soil water redistribution
    if SPARTREE == 1
        TEtree0_imp = TEtree * Ctree * (4 * r_tree) * Cimp
        TEtree0_bare = TEtree * Ctree * (4 * r_tree) * Cbare
        TEtree0_veg = TEtree * Ctree * (4 * r_tree) * Cveg
    else
        if (4 * r_tree) <= (fveg + fbare)
            TEtree0_imp = zero(FT) * Ctree * Cimp
            TEtree0_bare = TEtree * (4 * r_tree) / (fveg + fbare) * Ctree * Cbare
            TEtree0_veg = TEtree * (4 * r_tree) / (fveg + fbare) * Ctree * Cveg
        else
            TEtree0_imp = ((4 * r_tree) - (fveg + fbare)) * TEtree / fimp * Ctree * Cimp
            TEtree0_bare = TEtree * Ctree * Cbare
            TEtree0_veg = TEtree * Ctree * Cveg
        end
    end

    TEgveg = TEgveg * Cveg
    Egbare_soil = Egbare_soil * Cbare
    Egveg_soil = Egveg_soil * Cveg

    # Impervious soil column (from layer 3)
    V_gimp1, O_gimp1, _, Lk_gimp1, _, _, Exwat_H_gimp1, Exwat_L_gimp1, Rd_gimp1, TEgveg_imp1, TEtree_imp1, Egimp_soil1, _, WB_Soil_gimp1, Psi_Soil_gimp1, Kf_gimp1 = water_soil(
        Otm1_imp[3:end],
        f_inf_gimp,
        TEtree0_imp,
        zero(FT),
        zero(FT),
        zeros(FT, length(Qin_imp_tm1[3:end])),
        dth,
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
        [zero(FT)],
        Rrootl_H,
        [zero(FT)],
        PsiL50_H,
        [zero(FT)],
        PsiX50_H,
        Zs[3:end] .- Zs[3],
        row,
    )

    V_gimp1 = [NaN; NaN; V_gimp1]

    # Bare soil column
    V_gbare1, _, _, Lk_gbare1, _, _, _, _, Rd_gbare1, TEgveg_bare1, TEtree_bare1, Egbare_Soil1, _, WB_Soil_gbare1, _, _ = water_soil(
        Otm1_bare,
        f_inf_gbare,
        TEtree0_bare,
        zero(FT),
        Egbare_soil,
        zeros(FT, length(Qin_bare_tm1)),
        dth,
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
        [zero(FT)],
        Rrootl_H,
        [zero(FT)],
        PsiL50_H,
        [zero(FT)],
        PsiX50_H,
        Zs,
        row,
    )

    # Vegetated soil column
    V_gveg1, _, _, Lk_gveg1, _, _, _, _, Rd_gveg1, TEgveg_veg1, TEtree_veg1, Egveg_Soil1, _, WB_Soil_gveg1, _, _ = water_soil(
        Otm1_veg,
        f_inf_gveg,
        TEtree0_veg,
        TEgveg,
        Egveg_soil,
        zeros(FT, length(Qin_veg_tm1)),
        dth,
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
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
        Zs,
        row,
    )

    # Lateral water redistribution parameters
    T_SPAN = (zero(FT), dth)
    OPT_SM = (abstol=FT(0.05), maxstep=dth)
    dz = diff(Zs)

    # Get soil parameters
    _, _, _, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, _, _, _, _, _, _, _, _, _, _ = Soil.soil_parameters_total(
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

    # Initialize redistributed volumes
    V_gimp2 = zeros(FT, length(dz))
    V_gbare2 = zeros(FT, length(dz))
    V_gveg2 = zeros(FT, length(dz))

    # First two layers: exchange between bare and vegetated only
    for i in 1:2
        Vlat1 = [V_gbare1[i], V_gveg1[i]]

        sol = solve(
            ODEProblem(
                (V, p, t) -> Soil.soil_moistures_rich_comp_lat2(
                    V,
                    dz[i],
                    SPAR,
                    Ks_Zs[i],
                    Osat[i],
                    Ohy[i],
                    L[i],
                    Pe[i],
                    O33[i],
                    alpVG[i],
                    nVG[i],
                    Cbare,
                    Cveg,
                    fbare,
                    fveg,
                    Wcan,
                ),
                Vlat1,
                T_SPAN,
            ),
            Rosenbrock23(; autodiff=AutoFiniteDiff());
            abstol=OPT_SM.abstol,
            dtmax=OPT_SM.maxstep,
        )

        V_gimp2[i] = NaN
        V_gbare2[i] = sol.u[end][1]
        V_gveg2[i] = sol.u[end][2]
    end

    # Remaining layers: exchange between all three columns
    for i in 3:length(dz)
        Vlat1 = [V_gimp1[i], V_gbare1[i], V_gveg1[i]]

        sol = solve(
            ODEProblem(
                (V, p, t) -> Soil.soil_moistures_rich_comp_lat3(
                    V,
                    dz[i],
                    SPAR,
                    Ks_Zs[i],
                    Osat[i],
                    Ohy[i],
                    L[i],
                    Pe[i],
                    O33[i],
                    alpVG[i],
                    nVG[i],
                    Cimp,
                    Cbare,
                    Cveg,
                    fimp,
                    fbare,
                    fveg,
                    Wcan,
                ),
                Vlat1,
                T_SPAN,
            ),
            Rosenbrock23(; autodiff=AutoFiniteDiff());
            abstol=OPT_SM.abstol,
            dtmax=OPT_SM.maxstep,
        )

        V_gimp2[i] = sol.u[end][1]
        V_gbare2[i] = sol.u[end][2]
        V_gveg2[i] = sol.u[end][3]
    end

    # Back compute lateral fluxes
    Qin_imp = V_gimp2 .- V_gimp1            # [mm/dth]
    Qin_bare = V_gbare2 .- V_gbare1         # [mm/dth]
    Qin_veg = V_gveg2 .- V_gveg1           # [mm/dth]

    Qin_bare2imp = zeros(FT, length(dz))
    Qin_bare2veg = zeros(FT, length(dz))
    Qin_imp2bare = zeros(FT, length(dz))
    Qin_imp2veg = zeros(FT, length(dz))
    Qin_veg2imp = zeros(FT, length(dz))
    Qin_veg2bare = zeros(FT, length(dz))

    Qin_bare2imp[1:2] .= NaN
    Qin_imp2bare[1:2] .= NaN
    Qin_imp2veg[1:2] .= NaN
    Qin_veg2imp[1:2] .= NaN

    # Impervious
    V_gimp2, O_gimp2, OS_gimp2, Psi_Soil_gimp2, Psi_s_H_gimp2, Psi_s_L_gimp2, Exwat_H_gimp2, Exwat_L_gimp2, Kf_gimp2 = Soil.soil_moisture_conductivity_update(
        V_gimp2[3:end],
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
        Zs[3:end] .- Zs[3],
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    V_gimp2 = [NaN; NaN; V_gimp2]
    O_gimp2 = [NaN; NaN; O_gimp2]
    Exwat_H_gimp2 = [NaN; NaN; vec(Exwat_H_gimp2)]
    Exwat_L_gimp2 = [NaN; NaN; vec(Exwat_L_gimp2)]
    Psi_Soil_gimp2 = [NaN; NaN; Psi_Soil_gimp2]
    Kf_gimp2 = [NaN; NaN; Kf_gimp2]

    # Bare
    V_gbare2, O_gbare2, OS_gbare2, Psi_soil_gbare2, Psi_s_H_gbare2, Psi_s_L_gbare2, Exwat_H_gbare2, Exwat_L_gbare2, Kf_gbare2 = Soil.soil_moisture_conductivity_update(
        V_gbare2,
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
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    # Vegetated
    V_gveg2, O_gveg2, OS_gveg2, Psi_soil_gveg2, Psi_s_H_gveg2, Psi_s_L_gveg2, Exwat_H_gveg2, Exwat_L_gveg2, Kf_gveg2 = Soil.soil_moisture_conductivity_update(
        V_gveg2,
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
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    # Change in water volume
    Vtm1_imp = (Otm1_imp - Ohy) .* dz
    Vtm1_imp[1:2] .= NaN
    Vtm1_bare = (Otm1_bare - Ohy) .* dz
    Vtm1_veg = (Otm1_veg - Ohy) .* dz

    function nansum(x)
        return sum(filter(!isnan, x))
    end

    function nansum(arr, dims)
        return sum(x -> !isnan(x) * x, arr; dims=dims)
    end

    dV_dt_gimpTot = nansum(V_gimp2) - nansum(Vtm1_imp)
    dV_dt_gbareTot = nansum(V_gbare2) - nansum(Vtm1_bare)
    dV_dt_gvegTot = nansum(V_gveg2) - nansum(Vtm1_veg)

    # Surface water balance: Tree, impervious, bare, vegetated ground
    WBsurf_tree =
        Ctree * Rain - Etree_In * dth * 3600 * 1000 / row - q_tree_dwn - dIn_tree_dt # Is not used [mm/dth]

    WBsurf_imp =
        Cimp * Rain_ground + Cimp * Runon_tm1 - Egimp_Pond * dth * 3600 * 1000 / row -
        q_gimp_runoff - f_inf_gimp * dth - dIn_gimp_dt # [mm/dth]

    WBsurf_bare =
        Cbare * Rain_ground + Cbare * Runon_tm1 + Anthropogenic.Waterf_canyonBare -
        Egbare_Pond * dth * 3600 * 1000 / row - f_inf_gbare * dth - q_gbare_runoff -
        dIn_gbare_dt # [mm/dth]

    WBsurf_veg =
        Cveg * Rain_ground + Cveg * Runon_tm1 + Anthropogenic.Waterf_canyonVeg -
        (Egveg_In + Egveg_Pond) * dth * 3600 * 1000 / row - f_inf_gveg * dth -
        q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt # [mm/dth]

    # Soil water balance: impervious, bare, vegetated ground
    WBsoil_imp =
        f_inf_gimp * dth + nansum(Qin_imp) .-
        (TEgveg_imp1 .+ TEtree_imp1 .+ Egimp_soil1) * dth * 3600 * 1000 / row -
        Lk_gimp1 * dth - Rd_gimp1 - dV_dt_gimpTot # [mm/dth]

    WBsoil_bare =
        f_inf_gbare * dth + nansum(Qin_bare) .-
        (TEgveg_bare1 .+ TEtree_bare1 + Egbare_Soil1) * dth * 3600 * 1000 / row -
        Lk_gbare1 * dth - Rd_gbare1 - dV_dt_gbareTot # [mm/dth]

    WBsoil_veg =
        f_inf_gveg * dth + nansum(Qin_veg) .-
        (TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1) * dth * 3600 * 1000 / row -
        Lk_gveg1 * dth - Rd_gveg1 - dV_dt_gvegTot # [mm/dth]

    # Ground water balance: impervious, bare, vegetated ground
    WBimp_tot =
        Cimp * Rain_ground + Cimp * Runon_tm1 + nansum(Qin_imp) -
        (Egimp_Pond + TEgveg_imp1 + TEtree_imp1 + Egimp_soil1) * dth * 3600 * 1000 / row -
        q_gimp_runoff - dIn_gimp_dt - Lk_gimp1 * dth - Rd_gimp1 - dV_dt_gimpTot # [mm/dth]

    WBbare_tot =
        Cbare * Rain_ground +
        Cbare * Runon_tm1 +
        nansum(Qin_bare) +
        Anthropogenic.Waterf_canyonBare -
        (Egbare_Pond + TEgveg_bare1 + TEtree_bare1 + Egbare_Soil1) * dth * 3600 * 1000 /
        row - q_gbare_runoff - dIn_gbare_dt - Lk_gbare1 * dth - Rd_gbare1 - dV_dt_gbareTot # [mm/dth]

    WBveg_tot =
        Cveg * Rain_ground +
        Cveg * Runon_tm1 +
        nansum(Qin_veg) +
        Anthropogenic.Waterf_canyonVeg -
        (Egveg_In + Egveg_Pond + TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1) *
        dth *
        3600 *
        1000 / row - q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt - Lk_gveg1 * dth -
        Rd_gveg1 - dV_dt_gvegTot # [mm/dth]

    # Calculate total runoff including saturated excess
    Runoff =
        Per_runoff * (
            fimp * (q_gimp_runoff + Rd_gimp1) +
            fbare * (q_gbare_runoff + Rd_gbare1) +
            fveg * (q_gveg_runoff + Rd_gveg1)
        )  # [mm/dth]

    # Calculate total runon
    Runon =
        (1 - Per_runoff) * (
            fimp * (q_gimp_runoff + Rd_gimp1) +
            fbare * (q_gbare_runoff + Rd_gbare1) +
            fveg * (q_gveg_runoff + Rd_gveg1)
        )  # [mm/dth]

    # Calculate total evapotranspiration and leakage
    Etot =
        (
            fimp * Egimp_Pond +
            fbare * (Egbare_Pond + Egbare_soil) +
            fveg * (Egveg_Pond + Egveg_soil + TEgveg) +
            TEtree
        ) *
        3600 *
        1000 / row  # [mm/h]

    # Deep ground leakage and change in storage
    DeepGLk = fimp * Lk_gimp1 + fbare * Lk_gbare1 + fveg * Lk_gveg1  # [mm/h]

    StorageTot =
        fimp * (dIn_gimp_dt + dV_dt_gimpTot) +
        fbare * (dIn_gbare_dt + dV_dt_gbareTot) +
        fveg * (dIn_gveg_dt + dIn_gveg_pond_dt + dV_dt_gvegTot) +
        (4 * r_tree) * dIn_tree_dt  # [mm/dth]

    WBcanyon_flux =
        Rain +
        Runon_tm1 +
        fveg * Anthropogenic.Waterf_canyonVeg +
        fbare * Anthropogenic.Waterf_canyonBare - Etot - DeepGLk - Runoff - Runon -
        StorageTot  # [mm]

    WBtree_level =
        (4 * r_tree) *
        (Ctree * Rain - Etree_In * dth * 3600 * 1000 / row - q_tree_dwn - dIn_tree_dt)  # [mm/dth]

    WBground_level =
        fimp * (
            Rain_ground + Runon_tm1 - dIn_gimp_dt - Egimp_Pond * dth * 3600 * 1000 / row -
            q_gimp_runoff - Rd_gimp1 - f_inf_gimp * dth
        ) +
        fbare * (
            Rain_ground + Runon_tm1 + Anthropogenic.Waterf_canyonBare - dIn_gbare_dt -
            Egbare_Pond * dth * 3600 * 1000 / row - q_gbare_runoff - Rd_gbare1 -
            f_inf_gbare * dth
        ) +
        fveg * (
            Rain_ground + Runon_tm1 + Anthropogenic.Waterf_canyonVeg - dIn_gveg_dt -
            dIn_gveg_pond_dt - (Egveg_In + Egveg_Pond) * dth * 3600 * 1000 / row -
            q_gveg_runoff - Rd_gveg1 - f_inf_gveg * dth
        )  # [mm/dth]

    WBsoil_level =
        fimp * (
            f_inf_gimp * dth +
            3600 * 1000 / row * dth * (-TEgveg_imp1 - TEtree_imp1 - Egimp_soil1) -
            Lk_gimp1 * dth - dV_dt_gimpTot
        ) +
        fbare * (
            f_inf_gbare * dth +
            3600 * 1000 / row * dth * (-TEgveg_bare1 - TEtree_bare1 - Egbare_Soil1) -
            Lk_gbare1 * dth - dV_dt_gbareTot
        ) +
        fveg * (
            f_inf_gveg * dth +
            3600 * 1000 / row * dth * (-TEgveg_veg1 - TEtree_veg1 - Egveg_Soil1) -
            Lk_gveg1 * dth - dV_dt_gvegTot
        )  # [mm/dth]

    WBcanyon_level =
        Rain +
        Runon_tm1 +
        (4 * r_tree) * (-Etree_In * dth * 3600 * 1000 / row - dIn_tree_dt) +
        fimp * (
            -dIn_gimp_dt - Egimp_Pond * dth * 3600 * 1000 / row - q_gimp_runoff - Rd_gimp1 -
            3600 * 1000 * dth / row * (TEgveg_imp1 + TEtree_imp1 + Egimp_soil1) -
            Lk_gimp1 * dth - dV_dt_gimpTot
        ) +
        fbare * (
            Anthropogenic.Waterf_canyonBare - dIn_gbare_dt -
            Egbare_Pond * dth * 3600 * 1000 / row - q_gbare_runoff - Rd_gbare1 -
            3600 * 1000 / row * dth * (TEgveg_bare1 + TEtree_bare1 + Egbare_Soil1) -
            Lk_gbare1 * dth - dV_dt_gbareTot
        ) +
        fveg * (
            Anthropogenic.Waterf_canyonVeg - dIn_gveg_dt - dIn_gveg_pond_dt -
            (Egveg_In + Egveg_Pond) * dth * 3600 * 1000 / row - q_gveg_runoff - Rd_gveg1 -
            3600 * 1000 / row * dth * (TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1) -
            Lk_gveg1 * dth - dV_dt_gvegTot
        )  # [mm/dth]

    # Assign final variables
    V_gimp = V_gimp2
    O_gimp = O_gimp2
    V_gimp[isnan.(V_gimp)] .= 0
    O_gimp[isnan.(O_gimp)] .= 0
    OS_gimp = OS_gimp2
    Lk_gimp = Lk_gimp1
    Psi_s_H_gimp = Psi_s_H_gimp2
    Psi_s_L_gimp = Psi_s_L_gimp2
    Exwat_H_gimp = vec(Exwat_H_gimp2)
    Exwat_L_gimp = vec(Exwat_L_gimp2)
    Exwat_H_gimp[isnan.(Exwat_H_gimp)] .= 0
    Exwat_L_gimp[isnan.(Exwat_L_gimp)] .= 0
    Rd_gimp = Rd_gimp1
    TEgveg_imp = TEgveg_imp1
    Egimp_soil = Egimp_soil1
    dV_dt_gimp = dV_dt_gimpTot
    Psi_soil_gimp = Psi_Soil_gimp2
    Kf_gimp = Kf_gimp2

    # Update bare soil variables
    V_gbare = V_gbare2
    O_gbare = O_gbare2
    OS_gbare = OS_gbare2
    Lk_gbare = Lk_gbare1
    Psi_s_H_gbare = Psi_s_H_gbare2
    Psi_s_L_gbare = Psi_s_L_gbare2
    Exwat_H_gbare = Exwat_H_gbare2
    Exwat_L_gbare = Exwat_L_gbare2
    Rd_gbare = Rd_gbare1
    TEgveg_bare = TEgveg_bare1
    Egbare_Soil = Egbare_Soil1
    dV_dt_gbare = dV_dt_gbareTot
    Psi_soil_gbare = Psi_soil_gbare2
    Kf_gbare = Kf_gbare2

    # Update vegetated soil variables
    V_gveg = V_gveg2
    O_gveg = O_gveg2
    OS_gveg = OS_gveg2
    Lk_gveg = Lk_gveg1
    Psi_s_H_gveg = Psi_s_H_gveg2
    Psi_s_L_gveg = Psi_s_L_gveg2
    Exwat_H_gveg = Exwat_H_gveg2
    Exwat_L_gveg = Exwat_L_gveg2
    Rd_gveg = Rd_gveg1
    TEgveg_veg = TEgveg_veg1
    Egveg_Soil = Egveg_Soil1
    dV_dt_gveg = dV_dt_gvegTot
    Psi_soil_gveg = Psi_soil_gveg2
    Kf_gveg = Kf_gveg2

    # Update water balance variables
    WB_Soil_gimp = WB_Soil_gimp1
    WB_Soil_gbare = WB_Soil_gbare1
    WB_Soil_gveg = WB_Soil_gveg1

    # Update total values and rescale
    V = nansum(hcat(fimp * V_gimp, fbare * V_gbare, fveg * V_gveg), 2)
    O = nansum(hcat(fimp * O_gimp, fbare * O_gbare, fveg * O_gveg), 2)
    OS = nansum(hcat(fimp * OS_gimp, fbare * OS_gbare, fveg * OS_gveg), 2)
    Lk = nansum(hcat(fimp * Lk_gimp, fbare * Lk_gbare, fveg * Lk_gveg), 2)
    Rd = nansum(hcat(fimp * Rd_gimp, fbare * Rd_gbare, fveg * Rd_gveg), 2)
    dV_dt = nansum(hcat(fimp * dV_dt_gimp, fbare * dV_dt_gbare, fveg * dV_dt_gveg), 2)

    # Rescale first two layers
    V[1:2] ./= (fbare + fveg)
    O[1:2] ./= (fbare + fveg)
    V[isnan.(V)] .= 0
    O[isnan.(O)] .= 0

    # Assign vegetation parameters
    Psi_s_L = Psi_s_L_gveg
    Exwat_L = Exwat_L_gveg
    TEgveg_tot = TEgveg_veg

    Psi_s_H_tot = Psi_s_H_gveg
    Exwat_H = @. fimp * Exwat_H_gimp + fbare * Exwat_H_gbare + fveg * Exwat_H_gveg

    # Handle tree transpiration based on SPARTREE case
    if SPARTREE == 1  # Tree roots can access all water in the soil
        TEtree_imp = TEtree_imp1 / (4 * r_tree) * Cimp
        TEtree_bare = TEtree_bare1 / (4 * r_tree) * Cbare
        TEtree_veg = TEtree_veg1 / (4 * r_tree) * Cveg
    else  # Tree roots access based on crown size
        if (4 * r_tree) <= (fveg + fbare)
            TEtree_imp = 0 * Cimp
            TEtree_bare = TEtree_bare1 / (4 * r_tree) / (fveg + fbare) * Cbare
            TEtree_veg = TEtree_veg1 / (4 * r_tree) / (fveg + fbare) * Cveg
        else
            TEtree_imp = TEtree_imp1 / ((4 * r_tree) - (fveg + fbare)) * fimp * Cimp
            TEtree_bare = TEtree_bare1 * Cbare
            TEtree_veg = TEtree_veg1 * Cveg
        end
    end

    # Calculate total tree transpiration
    TEtree_tot = fimp * TEtree_imp + fbare * TEtree_bare + fveg * TEtree_veg

    # Compute energy balance terms
    EB_TEtree = TEtree_tot - TEtree
    EB_TEgveg = TEgveg_tot - TEgveg

    # Create water balance structs using NamedTuples
    WBIndv = (
        WB_In_tree=WB_In_tree,
        WB_In_gveg=WB_In_gveg,
        WB_In_gimp=WB_In_gimp,
        WB_In_gbare=WB_In_gbare,
        WB_Pond_gveg=WB_Pond_gveg,
        WB_Soil_gimp=WB_Soil_gimp,
        WB_Soil_gbare=WB_Soil_gbare,
        WB_Soil_gveg=WB_Soil_gveg,
    )

    WBTot = (
        WBsurf_tree=WBsurf_tree,
        WBsurf_imp=WBsurf_imp,
        WBsurf_bare=WBsurf_bare,
        WBsurf_veg=WBsurf_veg,
        WBsoil_imp=WBsoil_imp,
        WBsoil_bare=WBsoil_bare,
        WBsoil_veg=WBsoil_veg,
        WBimp_tot=WBimp_tot,
        WBbare_tot=WBbare_tot,
        WBveg_tot=WBveg_tot,
        WBcanyon_flux=WBcanyon_flux,
        WBtree_level=WBtree_level,
        WBground_level=WBground_level,
        WBsoil_level=WBsoil_level,
        WBcanyon_level=WBcanyon_level,
    )

    return (;
        q_tree_dwn,
        In_tree,
        dIn_tree_dt,
        q_gveg_dwn,
        In_gveg,
        dIn_gveg_dt,
        q_gimp_runoff,
        In_gimp,
        dIn_gimp_dt,
        f_inf_gimp,
        q_gbare_runoff,
        In_gbare,
        dIn_gbare_dt,
        f_inf_gbare,
        q_gveg_runoff,
        In_gveg_pond,
        dIn_gveg_pond_dt,
        f_inf_gveg,
        V_gimp,
        O_gimp,
        OS_gimp,
        Lk_gimp,
        Psi_s_H_gimp,
        Psi_s_L_gimp,
        Exwat_H_gimp,
        Exwat_L_gimp,
        Rd_gimp,
        TEgveg_imp,
        TEtree_imp,
        Egimp_soil,
        dV_dt_gimp,
        Psi_soil_gimp,
        Kf_gimp,
        V_gbare,
        O_gbare,
        OS_gbare,
        Lk_gbare,
        Psi_s_H_gbare,
        Psi_s_L_gbare,
        Exwat_H_gbare,
        Exwat_L_gbare,
        Rd_gbare,
        TEgveg_bare,
        TEtree_bare,
        Egbare_Soil,
        dV_dt_gbare,
        Psi_soil_gbare,
        Kf_gbare,
        V_gveg,
        O_gveg,
        OS_gveg,
        Lk_gveg,
        Psi_s_H_gveg,
        Psi_s_L_gveg,
        Exwat_H_gveg,
        Exwat_L_gveg,
        Rd_gveg,
        TEgveg_veg,
        TEtree_veg,
        Egveg_Soil,
        dV_dt_gveg,
        Psi_soil_gveg,
        Kf_gveg,
        Qin_imp,
        Qin_bare,
        Qin_veg,
        Qin_bare2imp,
        Qin_bare2veg,
        Qin_imp2bare,
        Qin_imp2veg,
        Qin_veg2imp,
        Qin_veg2bare,
        V,
        O,
        OS,
        Lk,
        Rd,
        dV_dt,
        Psi_s_L,
        Exwat_L,
        TEgveg_tot,
        Psi_s_H_tot,
        Exwat_H,
        TEtree_tot,
        EB_TEtree,
        EB_TEgveg,
        WBIndv,
        WBTot,
        Runoff,
        Runon_ittm,
        Etot,
        DeepGLk,
        StorageTot,
    )
end
