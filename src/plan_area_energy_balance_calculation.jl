
"""
    plan_area_energy_balance_calculation(
        results::Dict{String,Any},
        model::Model{FT},
        forcing::ModelComponents.ForcingInputSet{FT,1},
    ) where {FT<:AbstractFloat}

Calculates the plan area energy balance components for the urban area, canyon, and roof.
Analogous to `PlanAreaEnergyBalanceCalculation.m` in the original MATLAB code.

# Arguments
- `results`: Dictionary containing simulation results.
- `model`: The model structure.
- `forcing`: The forcing input data.

# Returns
- `EnergyFluxUrban`: Dictionary containing urban energy flux components.
- `EnergyFluxCan`: Dictionary containing canyon energy flux components.
- `EnergyFluxRoof`: Dictionary containing roof energy flux components.
"""
function plan_area_energy_balance_calculation(
    results::Dict{String,Any},
    model::Model{FT},
    forcing::ModelComponents.ForcingInputSet{FT,1},
    # figure::Bool,
    NN::Signed,
) where {FT<:AbstractFloat}

    # Incoming radiation
    Meteo = forcing.meteorological
    SWRin_atm =
        view(Meteo.SAB1_in, 1:NN) +
        view(Meteo.SAB2_in, 1:NN) +
        view(Meteo.SAD1_in, 1:NN) +
        view(Meteo.SAD2_in, 1:NN)
    SWRinDir_atm = view(Meteo.SAB1_in, 1:NN) + view(Meteo.SAB2_in, 1:NN)
    SWRinDiff_atm = view(Meteo.SAD1_in, 1:NN) + view(Meteo.SAD2_in, 1:NN)
    LWRin_atm = view(Meteo.LWR_in, 1:NN)

    # Unpack ViewFactors
    F_sg_T=results["ViewFactor"].F_sg_T
    F_sw_T=results["ViewFactor"].F_sw_T
    F_st_T=results["ViewFactor"].F_st_T

    # Normalize surface area
    urbangeometry = model.parameters.urbangeometry
    A_s = urbangeometry.wcanyon
    A_g = urbangeometry.wcanyon
    A_w = urbangeometry.hcanyon
    A_t = 2 * 2 * pi * urbangeometry.radius_tree

    fgveg = model.parameters.surfacefractions.ground.fveg
    fgbare = model.parameters.surfacefractions.ground.fbare
    fgimp = model.parameters.surfacefractions.ground.fimp

    wcanyon_norm = urbangeometry.wcanyon_norm
    wroof_norm = urbangeometry.wroof_norm
    Height_canyon = urbangeometry.hcanyon
    Width_canyon = urbangeometry.wcanyon
    Width_roof = urbangeometry.wroof

    BEM_on = true # TODO: Determine from model or results?

    # Rescale tree absorbed radiation
    SWRabsTree = results["SWRabsTree"] / FT(pi) # Technically 4*rad/(4*rad*pi)
    LWRabsTree = results["LWRabsTree"] / FT(pi)

    # Retrieve fluxes at this timestep
    # SWRabs = results["SWRabs_t"]
    # LWRabs = results["LWRabs_t"]

    # Roof specific fluxes (need to be in results)
    SWRinTotalRoof = results["SWRinTotalRoof"]
    SWRabsTotalRoof = results["SWRabsTotalRoof"]
    SWRoutTotalRoof = results["SWRoutTotalRoof"]
    LWRinTotalRoof = results["LWRinTotalRoof"]
    LWRabsTotalRoof = results["LWRabsTotalRoof"]
    LWRoutTotalRoof = results["LWRoutTotalRoof"]

    # Shortwave Canyon Components
    CanSWRin_SurfArea =
        A_g/A_g * (
            results["SWRinGroundVeg"] * fgveg +
            results["SWRinGroundBare"] * fgbare +
            results["SWRinGroundImp"] * fgimp
        ) +
        A_w/A_g * (results["SWRinWallSun"] + results["SWRinWallShade"]) +
        A_t/A_g * results["SWRinTree"]

    CanSWRabs_SurfArea =
        A_g/A_g * (
            results["SWRabsGroundVeg"] * fgveg +
            results["SWRabsGroundBare"] * fgbare +
            results["SWRabsGroundImp"] * fgimp
        ) +
        A_w/A_g * (results["SWRabsWallSun"] + results["SWRabsWallShade"]) +
        SWRabsTree * A_t/A_g

    CanSWRout_SurfArea =
        A_g/A_g * (
            results["SWRoutGroundVeg"] * fgveg +
            results["SWRoutGroundBare"] * fgbare +
            results["SWRoutGroundImp"] * fgimp
        ) +
        A_w/A_g * (results["SWRoutWallSun"] + results["SWRoutWallShade"]) +
        results["SWRoutTree"] * A_t/A_s

    CanSWRout_Ref_to_Atm =
        results["SWRoutGroundVeg"] * F_sg_T * fgveg +
        results["SWRoutGroundBare"] * F_sg_T * fgbare +
        results["SWRoutGroundImp"] * F_sg_T * fgimp +
        results["SWRoutWallSun"] * F_sw_T +
        results["SWRoutWallShade"] * F_sw_T +
        results["SWRoutTree"] * F_st_T

    CanSWR_EB_SurfArea = CanSWRin_SurfArea - CanSWRabs_SurfArea - CanSWRout_SurfArea
    CanSWR_EB_PlanAreaCanyon = SWRin_atm - CanSWRabs_SurfArea - CanSWRout_Ref_to_Atm

    # TODO: add missing EBwallSunSWR, EBwallShdSWR, SWREBinternal

    # Shortwave Urban Components
    UrbanSWRin_SurfArea = wcanyon_norm * CanSWRin_SurfArea + wroof_norm * SWRinTotalRoof
    UrbanSWRabs_SurfArea = wcanyon_norm * CanSWRabs_SurfArea + wroof_norm * SWRabsTotalRoof
    UrbanSWRout_SurfArea = wcanyon_norm * CanSWRout_SurfArea + wroof_norm * SWRoutTotalRoof
    UrbanSWRout_Ref_to_Atm =
        wcanyon_norm * CanSWRout_Ref_to_Atm + wroof_norm * SWRoutTotalRoof

    UrbanSWR_EB_SurfArea = UrbanSWRin_SurfArea - UrbanSWRabs_SurfArea - UrbanSWRout_SurfArea
    UrbanSWR_EB_PlanAreaUrban = SWRin_atm - UrbanSWRabs_SurfArea - UrbanSWRout_Ref_to_Atm

    # Canyon albedo calculation
    UrbanAlbedo = UrbanSWRout_Ref_to_Atm ./ SWRin_atm
    CanAlbedo = CanSWRout_Ref_to_Atm ./ SWRin_atm
    RoofAlbedo = model.parameters.optical.roof.albedo

    # Longwave Canyon Components
    CanLWRin_SurfArea =
        A_g/A_g * (
            results["LWRinGroundVeg"] * fgveg +
            results["LWRinGroundBare"] * fgbare +
            results["LWRinGroundImp"] * fgimp
        ) +
        A_w/A_g * (results["LWRinWallSun"] + results["LWRinWallShade"]) +
        results["LWRinTree"] * A_t/A_g

    CanLWRabs_SurfArea =
        A_g/A_g * (
            results["LWRabsGroundVeg"] * fgveg +
            results["LWRabsGroundBare"] * fgbare +
            results["LWRabsGroundImp"] * fgimp
        ) +
        A_w/A_g * (results["LWRabsWallSun"] + results["LWRabsWallShade"]) +
        LWRabsTree * A_t/A_g

    CanLWRout_SurfArea =
        A_g/A_s * (
            results["LWRoutGroundVeg"] * fgveg +
            results["LWRoutGroundBare"] * fgbare +
            results["LWRoutGroundImp"] * fgimp
        ) +
        A_w/A_s * (results["LWRoutWallSun"] + results["LWRoutWallShade"]) +
        results["LWRoutTree"] * A_t/A_s

    CanLWRout_Ref_to_Atm =
        results["LWRoutGroundVeg"] * F_sg_T * fgveg +
        results["LWRoutGroundBare"] * F_sg_T * fgbare +
        results["LWRoutGroundImp"] * F_sg_T * fgimp +
        results["LWRoutWallSun"] * F_sw_T +
        results["LWRoutWallShade"] * F_sw_T +
        results["LWRoutTree"] * F_st_T

    CanLWR_EB_SurfArea = CanLWRin_SurfArea - CanLWRabs_SurfArea - CanLWRout_SurfArea
    CanLWR_EB_PlanAreaCanyon = LWRin_atm - CanLWRabs_SurfArea - CanLWRout_Ref_to_Atm

    # Longwave Urban Components
    UrbanLWRin_SurfArea = wcanyon_norm * CanLWRin_SurfArea + wroof_norm * LWRinTotalRoof
    UrbanLWRabs_SurfArea = wcanyon_norm * CanLWRabs_SurfArea + wroof_norm * LWRabsTotalRoof
    UrbanLWRout_SurfArea = wcanyon_norm * CanLWRout_SurfArea + wroof_norm * LWRoutTotalRoof
    UrbanLWRout_Ref_to_Atm =
        wcanyon_norm * CanLWRout_Ref_to_Atm + wroof_norm * LWRoutTotalRoof

    UrbanLWR_EB_SurfArea = UrbanLWRin_SurfArea - UrbanLWRabs_SurfArea - UrbanLWRout_SurfArea
    UrbanLWR_EB_PlanAreaUrban = LWRin_atm - UrbanLWRabs_SurfArea - UrbanLWRout_Ref_to_Atm

    SWRabs_Urban = results["SWRabsTotalUrban"]
    LWRabs_Urban = results["LWRabsTotalUrban"]
    LE_Urban = results["LEfluxUrban"]
    H_Urban = results["HfluxUrban"]

    # LE Effect
    LEfluxRoof = results["LEfluxRoof"]
    LEfluxCanyon = results["LEfluxCanyon"] # Assuming added to results

    HfluxRoof = results["HfluxRoof"]
    HfluxCanyon = results["HfluxCanyon"]

    # Gflux
    G1Ground = results["G1Ground"] # Canyon ground

    # Retrieve struct/namedtuples
    Gfloor = results["GbuildIntGfloor"]     # Building floor
    G1WallSun = results["G1WallSun"]
    G1WallShade = results["G1WallShade"]
    G1Roof = results["G1Roof"]

    GdSinternalMass = results["GbuildIntdSinternalMass"]

    dsWallSun = results["dsWallSun"]
    dsWallShade = results["dsWallShade"]
    dsRoof = results["dsRoof"]

    if BEM_on
        Gground_Urban = wcanyon_norm * G1Ground + wroof_norm * Gfloor
        dSdt_buildEnv =
            wcanyon_norm * (A_w/A_g * (dsWallSun + dsWallShade)) +
            wroof_norm * dsRoof +
            Height_canyon / (Width_canyon + Width_roof) * GdSinternalMass
        G1Building = zero(FT)
    else
        Gground_Urban = wcanyon_norm * G1Ground
        G1Building =
            wcanyon_norm * (A_w/A_g * (G1WallSun + G1WallShade)) + wroof_norm * G1Roof
        dSdt_buildEnv = zero(FT)
    end

    dS_H_air = results["dS_H_air"]
    dS_LE_air = results["dS_LE_air"]

    # Roof internal?
    dSH_air_build = results["HbuildIntdSH_air"]
    dSLE_air_build = results["LEbuildIntdSLE_air"]

    dSdt_Air =
        wcanyon_norm * (dS_H_air + dS_LE_air) +
        wroof_norm * (dSH_air_build + dSLE_air_build)

    # Anthropogenic
    Qanth_Canyon = view(forcing.anthropogenic.Qf_canyon, 1:NN)
    Qanth_Roof = view(forcing.anthropogenic.Qf_roof, 1:NN)
    BEM_TotAnthInput_URB = results["WasteHeatTotAnthInput_URB"]

    Qanth = Qanth_Canyon * wcanyon_norm + Qanth_Roof * wroof_norm + BEM_TotAnthInput_URB

    BEM_WaterFromAC_Can = results["WasteHeatWaterFromAC_Can"]
    QanthACcondensation = BEM_WaterFromAC_Can * wcanyon_norm

    if BEM_on
        EBtot =
            SWRabs_Urban + LWRabs_Urban + Qanth + QanthACcondensation - LE_Urban - H_Urban -
            Gground_Urban - dSdt_buildEnv - dSdt_Air
    else
        EBtot =
            SWRabs_Urban + LWRabs_Urban + Qanth - LE_Urban - H_Urban - Gground_Urban -
            G1Building - dSdt_Air
    end

    EnergyFluxUrban = DataFrame(;
        SWRin_PlanArea=SWRin_atm,
        SWRin_SurfArea=UrbanSWRin_SurfArea,
        SWRabs_SurfArea=UrbanSWRabs_SurfArea,
        SWRout_SurfArea=UrbanSWRout_SurfArea,
        SWRout_Ref_to_Atm=UrbanSWRout_Ref_to_Atm,
        SWREB_SurfArea=UrbanSWR_EB_SurfArea,
        SWREB_PlanAreaUrban=UrbanSWR_EB_PlanAreaUrban,
        LWRin_PlanArea=LWRin_atm,
        LWRin_SurfArea=UrbanLWRin_SurfArea,
        LWRabs_SurfArea=UrbanLWRabs_SurfArea,
        LWRout_SurfArea=UrbanLWRout_SurfArea,
        LWRout_Ref_to_Atm=UrbanLWRout_Ref_to_Atm,
        LWREB_SurfArea=UrbanLWR_EB_SurfArea,
        LWREB_PlanAreaUrban=UrbanLWR_EB_PlanAreaUrban,
        UrbanAlbedo=UrbanAlbedo,
        SWRabs=SWRabs_Urban,
        LWRabs=LWRabs_Urban,
        LEflux=LE_Urban,
        Hflux=H_Urban,
        GfluxGround=Gground_Urban,
        G1Building=G1Building,
        dSdtBuildEnv=dSdt_buildEnv,
        dSdt_Air=dSdt_Air,
        Qanth=Qanth,
        QanthACcondensation=QanthACcondensation,
        EB=EBtot,
    )

    # Recalculate Transmitted components
    ParWindows = model.parameters.building_energy.windows
    GlazingRatio = ParWindows.GlazingRatio
    SWRabsWallSunTransmitted = GlazingRatio * results["SWRabsWallSun"]
    SWRabsWallShadeTransmitted = GlazingRatio * results["SWRabsWallShade"]
    BEM_SensibleFromVent_Can = results["WasteHeatSensibleFromVent_Can"]
    BEM_SensibleFromAC_Can = results["WasteHeatSensibleFromAC_Can"]
    BEM_SensibleFromHeat_Can = results["WasteHeatSensibleFromHeat_Can"]
    BEM_LatentFromVent_Can = results["WasteHeatLatentFromVent_Can"]
    BEM_LatentFromAC_Can = results["WasteHeatLatentFromAC_Can"]
    BEM_LatentFromHeat_Can = results["WasteHeatLatentFromHeat_Can"]

    Qanth =
        Qanth_Canyon +
        BEM_SensibleFromVent_Can +
        BEM_SensibleFromAC_Can +
        BEM_SensibleFromHeat_Can +
        BEM_LatentFromVent_Can +
        BEM_LatentFromAC_Can +
        BEM_LatentFromHeat_Can

    EnergyFluxCan = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        SWRin_SurfArea=CanSWRin_SurfArea,
        SWRabs_SurfArea=CanSWRabs_SurfArea,
        SWRout_SurfArea=CanSWRout_SurfArea,
        SWRout_Ref_to_Atm=CanSWRout_Ref_to_Atm,
        SWREB_SurfArea=CanSWR_EB_SurfArea,
        SWREB_PlanAreaCanyon=CanSWR_EB_PlanAreaCanyon,
        LWRin_SurfArea=CanLWRin_SurfArea,
        LWRabs_SurfArea=CanLWRabs_SurfArea,
        LWRout_SurfArea=CanLWRout_SurfArea,
        LWRout_Ref_to_Atm=CanLWRout_Ref_to_Atm,
        LWREB_SurfArea=CanLWR_EB_SurfArea,
        LWREB_PlanAreaCanyon=CanLWR_EB_PlanAreaCanyon,
        CanAlbedo=CanAlbedo,
        SWRabs=(
            results["SWRabsTotalCanyon"] -
            A_w/A_g * (SWRabsWallSunTransmitted + SWRabsWallShadeTransmitted)
        ),
        LWRabs=results["LWRabsTotalCanyon"],
        LEflux=LEfluxCanyon,
        Hflux=HfluxCanyon,
        Gflux=results["G1Canyon"],
        dSdt_Air=dS_H_air + dS_LE_air,
        Qanth=Qanth,
        QanthACcondensation=BEM_WaterFromAC_Can,
    )

    EnergyFluxCan[:, :EB] =
        EnergyFluxCan[:, "SWRabs"] +
        EnergyFluxCan[:, "LWRabs"] +
        EnergyFluxCan[:, "Qanth"] +
        EnergyFluxCan[:, "QanthACcondensation"] - EnergyFluxCan[:, "LEflux"] -
        EnergyFluxCan[:, "Hflux"] - EnergyFluxCan[:, "Gflux"] - EnergyFluxCan[:, "dSdt_Air"]

    # Fill EnergyFluxRoof
    EnergyFluxRoof = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        SWRin=SWRinTotalRoof,
        SWRout=SWRoutTotalRoof,
        LWRin=LWRinTotalRoof,
        LWRout=LWRoutTotalRoof,
        RoofAlbedo=RoofAlbedo,
        SWRabs=SWRabsTotalRoof,
        LWRabs=LWRabsTotalRoof,
        LEflux=LEfluxRoof,
        Hflux=HfluxRoof,
        Gflux=G1Roof,
        Qanth=Qanth_Roof,
        EB=SWRabsTotalRoof + LWRabsTotalRoof - LEfluxRoof - HfluxRoof - G1Roof + Qanth_Roof,
        SWREB=SWRinTotalRoof - SWRabsTotalRoof - SWRoutTotalRoof,
        LWREB=LWRinTotalRoof - LWRabsTotalRoof - LWRoutTotalRoof,
    )

    # if figure
    TTUrban = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        SWRin=EnergyFluxUrban[:, "SWRin_PlanArea"],
        SWRabs=EnergyFluxUrban[:, "SWRabs_SurfArea"],
        SWRout=EnergyFluxUrban[:, "SWRout_Ref_to_Atm"],
        LWRin=EnergyFluxUrban[:, "LWRin_PlanArea"],
        LWRabs=EnergyFluxUrban[:, "LWRabs_SurfArea"],
        LWRout=EnergyFluxUrban[:, "LWRout_Ref_to_Atm"],
        LE=EnergyFluxUrban[:, "LEflux"],
        H=EnergyFluxUrban[:, "Hflux"],
        G=EnergyFluxUrban[:, "GfluxGround"] +
          EnergyFluxUrban[:, "G1Building"] +
          EnergyFluxUrban[:, "dSdtBuildEnv"],
        Qanth=EnergyFluxUrban[:, "Qanth"],
        QantACcondensation=EnergyFluxUrban[:, "QanthACcondensation"],
        EB=EnergyFluxUrban[:, "EB"],
        Albedo=EnergyFluxUrban[:, "UrbanAlbedo"],
        BownRatio=EnergyFluxUrban[:, "Hflux"] ./ EnergyFluxUrban[:, "LEflux"],
    )

    TTUrbanDiurnal = @chain TTUrban begin
        select(Not(:Month))
        groupby(:Hour)
        combine(Not(:Hour) .=> NaNMath.mean, renamecols=false)
    end

    TTUrbanDiurnalMedian = @chain TTUrban begin
        select(Not(:Month))
        groupby(:Hour)
        combine(Not(:Hour) .=> NaNMath.median, renamecols=false)
    end

    TTUrbanSeasonal = @chain TTUrban begin
        select(Not(:Hour))
        groupby(:Month)
        combine(Not(:Month) .=> NaNMath.mean, renamecols=false)
    end

    TTUrbanSeasonalMedian = @chain TTUrban begin
        select(Not(:Hour))
        groupby(:Month)
        combine(Not(:Month) .=> NaNMath.median, renamecols=false)
    end

    fig1, fig2, fig3 = energy_balance_plots(
        TTUrban,
        TTUrbanDiurnal,
        TTUrbanDiurnalMedian,
        TTUrbanSeasonal,
        TTUrbanSeasonalMedian,
        view(forcing.datetime, 1:NN),
    )
    # end

    return EnergyFluxUrban, EnergyFluxCan, EnergyFluxRoof, fig1, fig2, fig3
end

function energy_budget_plots(
    EB_series::AbstractVector,
    Date::AbstractVector,
    TTUrbanDiurnal::DataFrame,
    TTUrbanSeasonal::DataFrame,
)
    # Time series
    p1 = plot(
        Date,
        EB_series;
        label="EB",
        color=:black,
        xlabel="time",
        ylabel="EB W/m²",
        title="Time series",
    )

    # Diurnal
    p2 = plot(
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.EB;
        label="EB",
        color=:black,
        linewidth=1.5,
        xlims=(0, 23),
        xlabel="hour",
        ylabel="EB W/m²",
        title="Diurnal",
    )

    # Seasonal
    p3 = plot(
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.EB;
        label="EB",
        color=:black,
        linewidth=1.5,
        xlims=(1, 12),
        xlabel="Month",
        ylabel="EB W/m²",
        title="Seasonal",
    )

    fig = plot(p1, p2, p3; layout=(1, 3), size=(1000, 500), title="Energy budget closure")
    return fig
end

function albedo_bowen_plots(
    TTUrbanDiurnal::DataFrame,
    TTUrbanDiurnalMedian::DataFrame,
    TTUrbanSeasonal::DataFrame,
    TTUrbanSeasonalMedian::DataFrame,
)
    # Diurnal Albedo
    p1 = plot(
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.Albedo;
        label="albedo",
        color=:black,
        linewidth=1.5,
        xlims=(0, 23),
        xlabel="hour",
        ylabel="Albedo (-)",
        title="Albedo = SWR_out/SWR_in\nDiurnal",
    )

    # Diurnal Bowen Ratio (Median)
    p2 = plot(
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnalMedian.BownRatio;
        label="Bowen ratio",
        color=:black,
        linewidth=1.5,
        xlims=(0, 23),
        xlabel="hour",
        ylabel="Bowen Ratio (-)",
        title="Bowen ratio = H/LE, median\nDiurnal",
    )

    # Seasonal Albedo
    p3 = plot(
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.Albedo;
        label="albedo",
        color=:black,
        linewidth=1.5,
        xlims=(1, 12),
        xlabel="month",
        ylabel="Albedo (-)",
        title="Seasonal",
    )

    # Seasonal Bowen Ratio (Median)
    p4 = plot(
        TTUrbanSeasonal.Month,
        TTUrbanSeasonalMedian.BownRatio;
        label="Bowen ratio",
        color=:black,
        linewidth=1.5,
        xlims=(1, 12),
        xlabel="month",
        ylabel="Bowen Ratio (-)",
        title="Seasonal",
    )

    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(800, 600))
    return fig
end

function energy_fluxes_plots(TTUrbanDiurnal::DataFrame, TTUrbanSeasonal::DataFrame)

    # Diurnal SWR
    p1 = plot(
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.SWRin;
        label="SWR_in",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p1,
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.SWRabs;
        label="SWR_abs",
        color=:red,
        linewidth=1.5,
    )
    plot!(
        p1,
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.SWRout;
        label="SWR_out",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p1;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="SWR W/m²",
        title="SWR\nDiurnal",
        legend=:bottom,
    )

    # Diurnal LWR
    p2 = plot(
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.LWRin;
        label="LWR_in",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p2,
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.LWRabs;
        label="LWR_abs",
        color=:red,
        linewidth=1.5,
    )
    plot!(
        p2,
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.LWRout;
        label="LWR_out",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p2;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="LWR W/m²",
        title="LWR\nDiurnal",
        legend=:bottom,
    )

    # Diurnal Energy Fluxes
    p3 = plot(
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.SWRabs;
        label="SWR_abs",
        color=:black,
        linewidth=1.5,
    )
    plot!(p3, TTUrbanDiurnal.Hour, TTUrbanDiurnal.H; label="H", color=:red, linewidth=1.5)
    plot!(
        p3, TTUrbanDiurnal.Hour, TTUrbanDiurnal.LE; label="LE", color=:green, linewidth=1.5
    )
    plot!(p3, TTUrbanDiurnal.Hour, TTUrbanDiurnal.G; label="G", color=:blue, linewidth=1.5)
    plot!(
        p3,
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.LWRabs;
        label="LWR_abs",
        color=:magenta,
        linewidth=1.5,
    )
    plot!(
        p3,
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.Qanth;
        label="Q_f",
        color=:black,
        linestyle=:dot,
        linewidth=1.5,
    )
    plot!(
        p3,
        TTUrbanDiurnal.Hour,
        TTUrbanDiurnal.QantACcondensation;
        label="Q_f,AC,cond",
        color=:black,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(
        p3;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="EB W/m²",
        title="Energy fluxes\nDiurnal",
        legend=:bottom,
        legend_columns=2,
    )

    # Seasonal SWR
    p4 = plot(
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.SWRin;
        label="SWR_in",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p4,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.SWRabs;
        label="SWR_abs",
        color=:red,
        linewidth=1.5,
    )
    plot!(
        p4,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.SWRout;
        label="SWR_out",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p4; xlims=(1, 12), xlabel="month", ylabel="SWR W/m²", title="Seasonal", legend=false
    )

    # Seasonal LWR
    p5 = plot(
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.LWRin;
        label="LWR_in",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p5,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.LWRabs;
        label="LWR_abs",
        color=:red,
        linewidth=1.5,
    )
    plot!(
        p5,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.LWRout;
        label="LWR_out",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p5; xlims=(1, 12), xlabel="month", ylabel="LWR W/m²", title="Seasonal", legend=false
    )

    # Seasonal Energy Fluxes
    p6 = plot(
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.SWRabs;
        label="SWR_abs",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p6, TTUrbanSeasonal.Month, TTUrbanSeasonal.H; label="H", color=:red, linewidth=1.5
    )
    plot!(
        p6,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.LE;
        label="LE",
        color=:green,
        linewidth=1.5,
    )
    plot!(
        p6, TTUrbanSeasonal.Month, TTUrbanSeasonal.G; label="G", color=:blue, linewidth=1.5
    )
    plot!(
        p6,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.LWRabs;
        label="LWR_abs",
        color=:magenta,
        linewidth=1.5,
    )
    plot!(
        p6,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.Qanth;
        label="Q_f",
        color=:black,
        linestyle=:dot,
        linewidth=1.5,
    )
    plot!(
        p6,
        TTUrbanSeasonal.Month,
        TTUrbanSeasonal.QantACcondensation;
        label="Q_f,AC,cond",
        color=:black,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(
        p6; xlims=(1, 12), xlabel="month", ylabel="EB W/m²", title="Seasonal", legend=false
    )

    fig = plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(1000, 800))
    return fig
end

function energy_balance_plots(
    TTUrban::DataFrame,
    TTUrbanDiurnal::DataFrame,
    TTUrbanDiurnalMedian::DataFrame,
    TTUrbanSeasonal::DataFrame,
    TTUrbanSeasonalMedian::DataFrame,
    Date::AbstractVector,
)
    fig1 = energy_budget_plots(TTUrban.EB, Date, TTUrbanDiurnal, TTUrbanSeasonal)
    fig2 = albedo_bowen_plots(
        TTUrbanDiurnal, TTUrbanDiurnalMedian, TTUrbanSeasonal, TTUrbanSeasonalMedian
    )
    fig3 = energy_fluxes_plots(TTUrbanDiurnal, TTUrbanSeasonal)

    return fig1, fig2, fig3
end
