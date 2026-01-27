"""
    water_balance_components(
        results::Dict{String,Any},
        model::Model{FT},
        forcing::ModelComponents.ForcingInputSet{FT,1},
        NN::Signed,
    ) where {FT<:AbstractFloat}

Calculates the water balance components for the urban area, canyon, roof, and building.
Analogous to `WaterBalanceComponents.m` in the original MATLAB code.

# Arguments
- `results`: Dictionary containing simulation results.
- `model`: The model structure.
- `forcing`: The forcing input data.
- `NN`: Number of timesteps.

# Returns
- `WaterFluxUrban`: DataFrame containing urban water flux components.
- `WaterFluxCan`: DataFrame containing canyon water flux components.
- `WaterFluxRoof`: DataFrame containing roof water flux components.
- `WaterFluxBuild`: DataFrame containing building water flux components.
"""
function water_balance_components(
    results::Dict{String,Any},
    model::Model{FT},
    forcing::ModelComponents.ForcingInputSet{FT,1},
    NN::Signed,
) where {FT<:AbstractFloat}

    # Calculation parameters
    Meteo = forcing.meteorological
    Tatm = view(Meteo.Tatm, 1:NN)
    L_heat = @. 1000.0 * (2501.3 - 2.361 * (Tatm - 273.15))  # Latent heat vaporization/condensation [J/kg]

    dts = 3600.0  # time step of calculation [s]
    dth = 1.0

    # Retrieve geometry and fractions
    urbangeometry = model.parameters.urbangeometry
    wcanyon_norm = urbangeometry.wcanyon_norm
    wroof_norm = urbangeometry.wroof_norm
    Width_canyon = urbangeometry.wcanyon
    Width_roof = urbangeometry.wroof
    radius_tree = urbangeometry.radius_tree

    frveg = model.parameters.surfacefractions.roof.fveg
    frimp = model.parameters.surfacefractions.roof.fimp

    fgveg = model.parameters.surfacefractions.ground.fveg
    fgbare = model.parameters.surfacefractions.ground.fbare
    fgimp = model.parameters.surfacefractions.ground.fimp

    # Precipitation
    RainRoof = view(Meteo.Rain, 1:NN)
    RainCan = view(Meteo.Rain, 1:NN)
    RainUrb = view(Meteo.Rain, 1:NN)

    # Irrigation at the surface
    IrrSurfRoof = frveg .* view(forcing.anthropogenic.Qf_roof, 1:NN)
    IrrSurfCan =
        fgveg .* view(forcing.anthropogenic.Waterf_canyonVeg, 1:NN) .+
        fgbare .* view(forcing.anthropogenic.Waterf_canyonBare, 1:NN)
    IrrSurfUrb = wcanyon_norm .* IrrSurfCan .+ wroof_norm .* IrrSurfRoof

    # Runoff leaving the system
    RunoffRoof = results["RunoffRoofTot"]
    RunoffCan = results["RunoffGroundTot"]
    RunoffUrb = results["RunoffUrban"]

    # Leakage at bottom of soil column
    LeakageRoof = dth .* results["LkRoof"]
    LeakageCan = dth .* results["LkGround"]
    LeakageUrb = dth .* results["LkUrban"]

    # Evapotranspiration, latent heat
    ETRoof = results["LEfluxRoof"] ./ L_heat .* dts
    ETCan = results["LEfluxCanyon"] ./ L_heat .* dts
    ETUrb = results["LEfluxUrban"] ./ L_heat .* dts

    # Change in latent heat stored in the air
    dS_ET_dtBuild = results["LEbuildIntdSLE_air"] ./ L_heat .* dts
    dS_ET_dtCan = results["dS_LE_air"] ./ L_heat .* dts
    dS_ET_dtUrb = wroof_norm .* dS_ET_dtBuild .+ wcanyon_norm .* dS_ET_dtCan

    # Urban evapotranspiration fluxes by source
    ETEvapoIntUrb =
        (
            wroof_norm .* (
                frimp .* results["LEfluxRoofImp"] .+
                frveg .* (results["LEfluxRoofVegInt"] .+ results["LEfluxRoofVegPond"])
            ) .+
            wcanyon_norm .* (
                fgimp .* results["LEfluxGroundImp"] .+
                fgbare .* results["LEfluxGroundBarePond"] .+
                fgveg .* (results["LEfluxGroundVegInt"] .+ results["LEfluxGroundVegPond"]) .+
                4.0 .* radius_tree .* results["LEfluxTreeInt"]
            )
        ) ./ L_heat .* dts

    ETEvapoSoilUrb =
        (
            wroof_norm .* frveg .* results["LEfluxRoofVegSoil"] .+
            wcanyon_norm .* (
                fgbare .* results["LEfluxGroundBareSoil"] .+
                fgveg .* results["LEfluxGroundVegSoil"]
            )
        ) ./ L_heat .* dts

    ETTranspUrb =
        (
            wroof_norm .* frveg .* results["LTEfluxRoofVeg"] .+
            wcanyon_norm .* (
                fgveg .* results["LTEfluxGroundVeg"] .+
                4.0 .* radius_tree .* results["LTEfluxTree"]
            )
        ) ./ L_heat .* dts

    ET_buildAnth =
        wroof_norm .* (results["LEbuildIntLEpeople"] .+ results["LEbuildIntLEequip"]) ./
        L_heat .* dts

    # Change in soil moisture
    dVdtRoof = frveg .* results["dVRoofSoilVeg_dt"]
    dVdtCan = results["dVGroundSoilTot_dt"]
    dVdtUrb = wcanyon_norm .* dVdtCan .+ wroof_norm .* dVdtRoof

    dVdtRoofCalc, dVdtCanCalc, dVdtUrbCalc = post_calculate_soil_moisture_change(
        results["OwaterInitial"],
        results,
        model.parameters.soil.roof,
        model.parameters.soil.ground,
        model.parameters.surfacefractions.roof,
        model.parameters.surfacefractions.ground,
        urbangeometry,
    )

    # Irrigation within soil (due to fixed soil moisture)
    IrrSoilRoof = dVdtRoofCalc .- dVdtRoof
    IrrSoilCan = dVdtCanCalc .- dVdtCan
    IrrSoilUrb = dVdtUrbCalc .- dVdtUrb

    # Change in intercepted water - On plant canopy
    dIdtPlantRoof = frveg .* results["dInt_dtRoofVegPlant"]
    dIdtPlantCan =
        fgveg .* results["dInt_dtGroundVegPlant"] .+
        4.0 .* radius_tree .* results["dInt_dtTree"]
    dIdtPlantUrb = wcanyon_norm .* dIdtPlantCan .+ wroof_norm .* dIdtPlantRoof

    # Change in intercepted water - On ground/surface
    dIdtGroundRoof =
        frveg .* results["dInt_dtRoofVegGround"] .+ frimp .* results["dInt_dtRoofImp"]
    dIdtGroundCan =
        fgveg .* results["dInt_dtGroundVegGround"] .+
        fgbare .* results["dInt_dtGroundBare"] .+ fgimp .* results["dInt_dtGroundImp"]
    dIdtGroundUrb = wcanyon_norm .* dIdtGroundCan .+ wroof_norm .* dIdtGroundRoof

    # Change in intercepted water - Due to runon
    dRun_dtRoof = vcat(results["RunonRoofTot"][1], diff(results["RunonRoofTot"]))
    dRun_dtCan = vcat(results["RunonGroundTot"][1], diff(results["RunonGroundTot"]))
    dRun_dtUrb = vcat(results["RunonUrban"][1], diff(results["RunonUrban"]))

    # Total change in intercepted water
    dIdtRoof = dIdtPlantRoof .+ dIdtGroundRoof .+ dRun_dtRoof
    dIdtCan = dIdtPlantCan .+ dIdtGroundCan .+ dRun_dtCan
    dIdtUrb = dIdtPlantUrb .+ dIdtGroundUrb .+ dRun_dtUrb

    # Surface water storage (SurfStor)
    IntRoof = results["IntRooftot"] .+ results["RunonRoofTot"]
    IntCan =
        fgimp .* results["IntGroundImp"] .+ fgbare .* results["IntGroundBare"] .+
        fgveg .* (results["IntGroundVegPlant"] .+ results["IntGroundVegGround"]) .+
        4.0 .* radius_tree .* results["IntTree"] .+ results["RunonGroundTot"]
    IntUrb = wcanyon_norm .* IntCan .+ wroof_norm .* IntRoof

    # Anthropogenic latent heat fluxes - Building interior sources
    ET_buildAnthRoof = zeros(FT, NN)
    ET_buildAnthCan = zeros(FT, NN)
    ET_buildAnthBuild =
        (results["LEbuildIntLEpeople"] .+ results["LEbuildIntLEequip"]) ./ L_heat .* dts
    ET_buildAnthUrb =
        wroof_norm .* (results["LEbuildIntLEpeople"] .+ results["LEbuildIntLEequip"]) ./
        L_heat .* dts

    # Latent heat removed from air due to ventilation
    ET_VentRoof = zeros(FT, NN)
    ET_VentCan = results["WasteHeatLatentFromVent_Can"] ./ L_heat .* dts
    ET_VentBuild = results["LEbuildIntLEvent"] ./ L_heat .* dts
    ET_VentUrb = zeros(FT, NN)

    # ET exchange due to AC
    ET_ACRoof = zeros(FT, NN)
    ET_ACCan = results["WasteHeatLatentFromAC_Can"] ./ L_heat .* dts
    ET_ACBuild =
        results["WasteHeatLatentFromAC_Can"] .* (Width_canyon ./ Width_roof) ./ L_heat .*
        dts
    ET_ACUrb = zeros(FT, NN)

    # ET exchange due to HVAC between indoor and outdoor air
    ET_HVACexchRoof = zeros(FT, NN)
    ET_HVACexchCan = ET_VentCan .+ ET_ACCan
    ET_HVACexchBuild = ET_VentBuild .- ET_ACBuild
    ET_HVACexchUrb = zeros(FT, NN)

    # Water removed from building interior due to condensation during AC
    ET_WasteWaterACBuild =
        results["WasteHeatWaterFromAC_Can"] .* (Width_canyon ./ Width_roof) ./ L_heat .* dts
    ET_WasteWaterACUrb =
        results["WasteHeatWaterFromAC_Can"] .* Width_canyon ./
        (Width_canyon .+ Width_roof) ./ L_heat .* dts

    # Water balance
    WBRoof =
        RainRoof .+ IrrSurfRoof .+ IrrSoilRoof .- RunoffRoof .- LeakageRoof .- ETRoof .-
        dVdtRoofCalc .- dIdtRoof
    WBCan =
        RainCan .+ IrrSurfCan .+ IrrSoilCan .- RunoffCan .- LeakageCan .- ETCan .-
        dVdtCanCalc .- dIdtCan .- dS_ET_dtCan
    WBBuild = ET_buildAnthBuild .+ ET_HVACexchBuild .- ET_WasteWaterACBuild .- dS_ET_dtBuild
    WBUrb =
        RainUrb .+ IrrSurfUrb .+ IrrSoilUrb .+ ET_buildAnthUrb .- RunoffUrb .- LeakageUrb .-
        ETUrb .- dVdtUrbCalc .- dIdtUrb .- dS_ET_dtUrb

    # Create DataFrames for output
    WaterFluxRoof = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        Rain=RainRoof,
        Runoff=RunoffRoof,
        Leakage=LeakageRoof,
        ET=ETRoof,
        dS_ET_dt=zeros(FT, NN),
        ET_HVACexch=ET_HVACexchRoof,
        AnthBuildInt=ET_buildAnthRoof,
        WasteWaterAC=zeros(FT, NN),
        dVdt=dVdtRoofCalc,
        dIdt=dIdtRoof,
        IrrSurf=IrrSurfRoof,
        IrrSoil=IrrSoilRoof,
        IrrTot=IrrSoilRoof .+ IrrSurfRoof,
        Int=IntRoof,
        WB=WBRoof,
    )

    WaterFluxCan = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        Rain=RainCan,
        Runoff=RunoffCan,
        Leakage=LeakageCan,
        ET=ETCan,
        dS_ET_dt=dS_ET_dtCan,
        ET_HVACexch=ET_HVACexchCan,
        AnthBuildInt=zeros(FT, NN),
        WasteWaterAC=zeros(FT, NN),
        dVdt=dVdtCanCalc,
        dIdt=dIdtCan,
        IrrSurf=IrrSurfCan,
        IrrSoil=IrrSoilCan,
        IrrTot=IrrSoilCan .+ IrrSurfCan,
        Int=IntCan,
        WB=WBCan,
    )

    WaterFluxBuild = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        dS_ET_dt=dS_ET_dtBuild,
        ET_HVACexch=ET_HVACexchBuild,
        AnthBuildInt=ET_buildAnthBuild,
        WasteWaterAC=ET_WasteWaterACBuild,
        WB=WBBuild,
    )

    WaterFluxUrban = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        Rain=RainUrb,
        Runoff=RunoffUrb,
        Leakage=LeakageUrb,
        ET=ETUrb,
        ETEvaporationFromSurface=ETEvapoIntUrb,
        ETEvaporationFromSoil=ETEvapoSoilUrb,
        ETTranspiration=ETTranspUrb,
        ETAnthSourceBuildInt=ET_buildAnth,
        dS_ET_dt=dS_ET_dtUrb,
        ET_HVACexch=ET_HVACexchUrb,
        AnthBuildInt=ET_buildAnthUrb,
        WasteWaterAC=ET_WasteWaterACUrb,
        dVdt=dVdtUrbCalc,
        dIdt=dIdtUrb,
        IrrSurf=IrrSurfUrb,
        IrrSoil=IrrSoilUrb,
        IrrTot=IrrSoilUrb .+ IrrSurfUrb,
        Int=IntUrb,
        WB=WBUrb,
    )

    # Compute diurnal and seasonal averages
    WaterFluxUrbanDiurnal = @chain WaterFluxUrban begin
        select(Not(:Month))
        groupby(:Hour)
        combine(Not(:Hour) .=> NaNMath.mean, renamecols=false)
    end

    WaterFluxUrbanSeasonal = @chain WaterFluxUrban begin
        select(Not(:Hour))
        groupby(:Month)
        combine(Not(:Month) .=> NaNMath.mean, renamecols=false)
    end

    # Generate plots
    fig1, fig2 = water_balance_plots(
        WaterFluxUrban,
        WaterFluxUrbanDiurnal,
        WaterFluxUrbanSeasonal,
        view(forcing.datetime, 1:NN),
    )

    return WaterFluxUrban, WaterFluxCan, WaterFluxRoof, WaterFluxBuild, fig1, fig2
end

function water_balance_plots(
    WaterFluxUrban::DataFrame,
    WaterFluxUrbanDiurnal::DataFrame,
    WaterFluxUrbanSeasonal::DataFrame,
    Date::AbstractVector,
)
    # Water budget WB
    p1 = plot(Date, WaterFluxUrban.WB; label="WB", color=:black, linewidth=1.5)
    plot!(p1; xlabel="Time", ylabel="WB (mm/time step)", title="Time series", grid=true)

    p2 = plot(
        WaterFluxUrbanDiurnal.Hour,
        WaterFluxUrbanDiurnal.WB;
        label="WB",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p2;
        xlims=(0, 23),
        xlabel="Hour",
        ylabel="WB (mm/time step)",
        title="Diurnal",
        grid=true,
    )

    p3 = plot(
        WaterFluxUrbanSeasonal.Month,
        WaterFluxUrbanSeasonal.WB;
        label="WB",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p3;
        xlims=(1, 12),
        xlabel="Month",
        ylabel="WB (mm/time step)",
        title="Seasonal",
        grid=true,
    )

    fig1 = plot(p1, p2, p3; layout=(1, 3), size=(1000, 300))

    # Evapotranspiration and interception
    p4 = plot(
        WaterFluxUrbanDiurnal.Hour,
        WaterFluxUrbanDiurnal.ET;
        label="ETₜₒₜ",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p4,
        WaterFluxUrbanDiurnal.Hour,
        WaterFluxUrbanDiurnal.ETEvaporationFromSurface;
        label="Eₛᵤᵣfₐcₑ",
        color=:red,
        linewidth=1.5,
    )
    plot!(
        p4,
        WaterFluxUrbanDiurnal.Hour,
        WaterFluxUrbanDiurnal.ETEvaporationFromSoil;
        label="Eₛₒᵢₗ",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p4,
        WaterFluxUrbanDiurnal.Hour,
        WaterFluxUrbanDiurnal.ETTranspiration;
        label="Tᵥₑgₑₜₐₜᵢₒₙ",
        color=:green,
        linewidth=1.5,
    )
    plot!(
        p4,
        WaterFluxUrbanDiurnal.Hour,
        WaterFluxUrbanDiurnal.AnthBuildInt;
        label="QLE,anth,building",
        color=:yellow,
        linewidth=1.5,
    )
    plot!(
        p4,
        WaterFluxUrbanDiurnal.Hour,
        -WaterFluxUrbanDiurnal.WasteWaterAC;
        label="QAC,waste water",
        color=:magenta,
        linewidth=1.5,
    )
    plot!(
        p4;
        xlims=(0, 23),
        xlabel="Hour",
        ylabel="ET (mm/time step)",
        title="Diurnal",
        grid=true,
    )

    p5 = plot(
        WaterFluxUrbanSeasonal.Month,
        WaterFluxUrbanSeasonal.ET;
        label="ETₜₒₜ",
        color=:black,
        linewidth=1.5,
    )
    plot!(
        p5,
        WaterFluxUrbanSeasonal.Month,
        WaterFluxUrbanSeasonal.ETEvaporationFromSurface;
        label="Eₛᵤᵣfₐcₑ",
        color=:red,
        linewidth=1.5,
    )
    plot!(
        p5,
        WaterFluxUrbanSeasonal.Month,
        WaterFluxUrbanSeasonal.ETEvaporationFromSoil;
        label="Eₛₒᵢₗ",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p5,
        WaterFluxUrbanSeasonal.Month,
        WaterFluxUrbanSeasonal.ETTranspiration;
        label="Tᵥₑgₑₜₐₜᵢₒₙ",
        color=:green,
        linewidth=1.5,
    )
    plot!(
        p5,
        WaterFluxUrbanSeasonal.Month,
        WaterFluxUrbanSeasonal.AnthBuildInt;
        label="QLE,anth,building",
        color=:yellow,
        linewidth=1.5,
    )
    plot!(
        p5,
        WaterFluxUrbanSeasonal.Month,
        -WaterFluxUrbanSeasonal.WasteWaterAC;
        label="QAC,waste water",
        color=:magenta,
        linewidth=1.5,
    )
    plot!(
        p5;
        xlims=(1, 12),
        xlabel="Month",
        ylabel="ET (mm/time step)",
        title="Seasonal",
        legend=:outertopright,
        grid=true,
    )

    fig2 = plot(p4, p5; layout=(1, 2), size=(1000, 400))

    return fig1, fig2
end
