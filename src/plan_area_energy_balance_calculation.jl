
"""
    plan_area_energy_balance_calculation(
        results::Dict{Symbol,Dict{Symbol,Array}},
        model::Model{FT},
        forcing::ModelComponents.ForcingInputSet{FT,1},
        view_factor::RayTracing.ViewFactor{FT},
        NN::Signed,
        BEM_on::Bool=true,
    ) where {FT<:AbstractFloat}

Calculates the plan area energy balance components for the urban area, canyon, and roof.
Analogous to `PlanAreaEnergyBalanceCalculation.m` in the original MATLAB code.

# Arguments
- `results`: Dictionary containing simulation results.
- `model`: The model structure.
- `forcing`: The forcing input data.
- `view_factor`: View factor structure for radiation calculations.
- `NN`: Number of time steps to consider.
- `BEM_on`: Boolean flag to indicate if Building Energy Model is active (default: true).

# Returns
- `EnergyFluxUrban`: Dictionary containing urban energy flux components.
- `EnergyFluxCan`: Dictionary containing canyon energy flux components.
- `EnergyFluxRoof`: Dictionary containing roof energy flux components.
- `fig1`, `fig2`, `fig3`: Figures to visualize the results.
"""
function plan_area_energy_balance_calculation(
    results::Dict{Symbol,Dict{Symbol,Array}},
    model::Model{FT},
    forcing::ModelComponents.ForcingInputSet{FT,1},
    view_factor::RayTracing.ViewFactor{FT},
    NN::Signed,
    BEM_on::Bool=true,
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
    F_sg_T=view_factor.F_sg_T
    F_sw_T=view_factor.F_sw_T
    F_st_T=view_factor.F_st_T

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

    # Rescale tree absorbed radiation
    SWRabsTree = results[:SWRabs][:Tree] / FT(pi) # Technically 4*rad/(4*rad*pi)
    LWRabsTree = results[:LWRabs][:Tree] / FT(pi)

    # Retrieve fluxes at this timestep
    # SWRabs = results["SWRabs_t"]
    # LWRabs = results["LWRabs_t"]

    # Roof specific fluxes (need to be in results)
    SWRinTotalRoof = results[:SWRin][:TotalRoof]
    SWRabsTotalRoof = results[:SWRabs][:TotalRoof]
    SWRoutTotalRoof = results[:SWRout][:TotalRoof]
    LWRinTotalRoof = results[:LWRin][:TotalRoof]
    LWRabsTotalRoof = results[:LWRabs][:TotalRoof]
    LWRoutTotalRoof = results[:LWRout][:TotalRoof]

    # Shortwave Canyon Components
    CanSWRin_SurfArea =
        A_g/A_g * (
            results[:SWRin][:GroundVeg] * fgveg +
            results[:SWRin][:GroundBare] * fgbare +
            results[:SWRin][:GroundImp] * fgimp
        ) +
        A_w/A_g * (results[:SWRin][:WallSun] + results[:SWRin][:WallShade]) +
        A_t/A_g * results[:SWRin][:Tree]

    CanSWRabs_SurfArea =
        A_g/A_g * (
            results[:SWRabs][:GroundVeg] * fgveg +
            results[:SWRabs][:GroundBare] * fgbare +
            results[:SWRabs][:GroundImp] * fgimp
        ) +
        A_w/A_g * (results[:SWRabs][:WallSun] + results[:SWRabs][:WallShade]) +
        SWRabsTree * A_t/A_g

    CanSWRout_SurfArea =
        A_g/A_g * (
            results[:SWRout][:GroundVeg] * fgveg +
            results[:SWRout][:GroundBare] * fgbare +
            results[:SWRout][:GroundImp] * fgimp
        ) +
        A_w/A_g * (results[:SWRout][:WallSun] + results[:SWRout][:WallShade]) +
        results[:SWRout][:Tree] * A_t/A_s

    CanSWRout_Ref_to_Atm =
        results[:SWRout][:GroundVeg] * F_sg_T * fgveg +
        results[:SWRout][:GroundBare] * F_sg_T * fgbare +
        results[:SWRout][:GroundImp] * F_sg_T * fgimp +
        results[:SWRout][:WallSun] * F_sw_T +
        results[:SWRout][:WallShade] * F_sw_T +
        results[:SWRout][:Tree] * F_st_T

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
            results[:LWRin][:GroundVeg] * fgveg +
            results[:LWRin][:GroundBare] * fgbare +
            results[:LWRin][:GroundImp] * fgimp
        ) +
        A_w/A_g * (results[:LWRin][:WallSun] + results[:LWRin][:WallShade]) +
        results[:LWRin][:Tree] * A_t/A_g

    CanLWRabs_SurfArea =
        A_g/A_g * (
            results[:LWRabs][:GroundVeg] * fgveg +
            results[:LWRabs][:GroundBare] * fgbare +
            results[:LWRabs][:GroundImp] * fgimp
        ) +
        A_w/A_g * (results[:LWRabs][:WallSun] + results[:LWRabs][:WallShade]) +
        LWRabsTree * A_t/A_g

    CanLWRout_SurfArea =
        A_g/A_s * (
            results[:LWRout][:GroundVeg] * fgveg +
            results[:LWRout][:GroundBare] * fgbare +
            results[:LWRout][:GroundImp] * fgimp
        ) +
        A_w/A_s * (results[:LWRout][:WallSun] + results[:LWRout][:WallShade]) +
        results[:LWRout][:Tree] * A_t/A_s

    CanLWRout_Ref_to_Atm =
        results[:LWRout][:GroundVeg] * F_sg_T * fgveg +
        results[:LWRout][:GroundBare] * F_sg_T * fgbare +
        results[:LWRout][:GroundImp] * F_sg_T * fgimp +
        results[:LWRout][:WallSun] * F_sw_T +
        results[:LWRout][:WallShade] * F_sw_T +
        results[:LWRout][:Tree] * F_st_T

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

    SWRabs_Urban = results[:SWRabs][:TotalUrban]
    LWRabs_Urban = results[:LWRabs][:TotalUrban]
    LE_Urban = results[:LEflux][:LEfluxUrban]
    H_Urban = results[:Hflux][:HfluxUrban]

    # LE Effect
    LEfluxRoof = results[:LEflux][:LEfluxRoof]
    LEfluxCanyon = results[:LEflux][:LEfluxCanyon] # Assuming added to results

    HfluxRoof = results[:Hflux][:HfluxRoof]
    HfluxCanyon = results[:Hflux][:HfluxCanyon]

    # Gflux
    G1Ground = results[:Gflux][:G1Ground] # Canyon ground

    # Retrieve struct/namedtuples
    Gfloor = results[:GbuildInt][:Gfloor]     # Building floor
    G1WallSun = results[:Gflux][:G1WallSun]
    G1WallShade = results[:Gflux][:G1WallShade]
    G1Roof = results[:Gflux][:G1Roof]

    GdSinternalMass = results[:GbuildInt][:dSinternalMass]

    dsWallSun = results[:dStorage][:dsWallSun]
    dsWallShade = results[:dStorage][:dsWallShade]
    dsRoof = results[:dStorage][:dsRoof]

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

    dS_H_air = results[:Hflux][:dS_H_air]
    dS_LE_air = results[:LEflux][:dS_LE_air]

    # Roof internal?
    dSH_air_build = results[:HbuildInt][:dSH_air]
    dSLE_air_build = results[:LEbuildInt][:dSLE_air]

    dSdt_Air =
        wcanyon_norm * (dS_H_air + dS_LE_air) +
        wroof_norm * (dSH_air_build + dSLE_air_build)

    # Anthropogenic
    Qanth_Canyon = view(forcing.anthropogenic.Qf_canyon, 1:NN)
    Qanth_Roof = view(forcing.anthropogenic.Qf_roof, 1:NN)
    BEM_TotAnthInput_URB = results[:BEMWasteHeat][:TotAnthInput_URB]

    Qanth = Qanth_Canyon * wcanyon_norm + Qanth_Roof * wroof_norm + BEM_TotAnthInput_URB

    BEM_WaterFromAC_Can = results[:BEMWasteHeat][:WaterFromAC_Can]
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
    SWRabsWallSunTransmitted = GlazingRatio * results[:SWRabs][:WallSun]
    SWRabsWallShadeTransmitted = GlazingRatio * results[:SWRabs][:WallShade]
    BEM_SensibleFromVent_Can = results[:BEMWasteHeat][:SensibleFromVent_Can]
    BEM_SensibleFromAC_Can = results[:BEMWasteHeat][:SensibleFromAC_Can]
    BEM_SensibleFromHeat_Can = results[:BEMWasteHeat][:SensibleFromHeat_Can]
    BEM_LatentFromVent_Can = results[:BEMWasteHeat][:LatentFromVent_Can]
    BEM_LatentFromAC_Can = results[:BEMWasteHeat][:LatentFromAC_Can]
    BEM_LatentFromHeat_Can = results[:BEMWasteHeat][:LatentFromHeat_Can]

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
            results[:SWRabs][:TotalCanyon] -
            A_w/A_g * (SWRabsWallSunTransmitted + SWRabsWallShadeTransmitted)
        ),
        LWRabs=results[:LWRabs][:TotalCanyon],
        LEflux=LEfluxCanyon,
        Hflux=HfluxCanyon,
        Gflux=results[:Gflux][:G1Canyon],
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
        QanthACcondensation=EnergyFluxUrban[:, "QanthACcondensation"],
        EB=EnergyFluxUrban[:, "EB"],
        Albedo=EnergyFluxUrban[:, "UrbanAlbedo"],
        BowenRatio=EnergyFluxUrban[:, "Hflux"] ./ EnergyFluxUrban[:, "LEflux"],
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
        TTUrbanDiurnalMedian.BowenRatio;
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
        TTUrbanSeasonalMedian.BowenRatio;
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
        TTUrbanDiurnal.QanthACcondensation;
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
        TTUrbanSeasonal.QanthACcondensation;
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
