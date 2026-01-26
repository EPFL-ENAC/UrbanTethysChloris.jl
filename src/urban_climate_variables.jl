"""
    urban_climate_variables(
        results::Dict{String,Any},
        model::Model{FT},
        forcing::ForcingInputSet{FT,1},
        NN::Signed,
    ) where {FT<:AbstractFloat}

Processes the simulation `results` to extract urban climate variables such as
air temperature, surface temperatures, relative humidity, mean radiant temperature,
and UTCI. It computes diurnal and seasonal averages and generates plots for
visualization.

# Arguments
- `results::Dict{String,Any}`: Dictionary containing simulation results.
- `model::Model{FT}`: The urban climate model used for the simulation.
- `forcing::ForcingInputSet{FT,1}`: Forcing input data used in the simulation.
- `NN::Signed`: Number of timesteps in the simulation.
"""
function urban_climate_variables(
    results::Dict{String,Any}, model::Model{FT}, forcing::ForcingInputSet{FT,1}, NN::Signed
) where {FT<:AbstractFloat}
    TTUrban = DataFrame(;
        Hour=Dates.hour.(view(forcing.datetime, 1:NN)),
        Month=Dates.month.(view(forcing.datetime, 1:NN)),
        T2m=tocelsius.(results["T2m"]),
        TCanyon=tocelsius.(results["TCanyon"]),
        RH2m=results["RH_T2m"] * 100,
        UTCI=results["UTCI"],
        Tmrt=results["Tmrt"],
        TRoofImp=tocelsius.(results["TRoofImp"]),
        TRoofVeg=tocelsius.(results["TRoofVeg"]),
        TGroundImp=tocelsius.(results["TGroundImp"]),
        TGroundBare=tocelsius.(results["TGroundBare"]),
        TGroundVeg=tocelsius.(results["TGroundVeg"]),
        TTree=tocelsius.(results["TTree"]),
        TWallSun=tocelsius.(results["TWallSun"]),
        TWallShade=tocelsius.(results["TWallShade"]),
        Tatm=tocelsius.(view(forcing.meteorological.Tatm, 1:NN)),
        RHatm=view(forcing.meteorological.rel_hum, 1:NN) * 100,
        Tbin=tocelsius.(results["Tbin"]),
        RHbin=results["RHbin"] * 100,
    )

    TTUrbanDiurnal = @chain TTUrban begin
        select(Not(:Month))
        groupby(:Hour)
        combine(Not(:Hour) .=> NaNMath.mean, renamecols=false)
    end

    TTUrbanSeasonal = @chain TTUrban begin
        select(Not(:Hour))
        groupby(:Month)
        combine(Not(:Month) .=> NaNMath.mean, renamecols=false)
    end

    fig1 = outdoor_thermal_comfort(TTUrbanDiurnal, TTUrbanSeasonal)
    fig2 = surface_temperatures(TTUrbanDiurnal, TTUrbanSeasonal, model)

    return TTUrban, TTUrbanDiurnal, TTUrbanSeasonal, fig1, fig2
end

const ZEROKELVIN = 273.15

function tocelsius(T::AbstractFloat)
    return T - ZEROKELVIN
end

function tokelvin(T::AbstractFloat)
    return T + ZEROKELVIN
end

function diurnal_UTCI(x::DataFrame)
    p = plot(x.Hour, x.UTCI; label="UTCI", color=:red, linewidth=1.5)
    plot!(
        p, x.Hour, (x.T2m); label=L"$T_{2m}$", color=:blue, linestyle=:dash, linewidth=1.5
    )
    plot!(
        p;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="T (°C)",
        title="Outdoor thermal comfort\nDiurnal",
        legend=:bottom,
        legend_columns=2,
        grid=true,
    )

    return p
end

function diurnal_Tmrt(x::DataFrame)
    p = plot(x.Hour, (x.Tmrt); label=L"$T_{mrt}$", color=:red, linewidth=1.5)
    plot!(
        p, x.Hour, (x.T2m); label=L"$T_{2m}$", color=:blue, linestyle=:dash, linewidth=1.5
    )
    plot!(
        p;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="T (°C)",
        title="Mean radiant temperature\nDiurnal",
        legend=:bottom,
        legend_columns=2,
        grid=true,
    )

    return p
end

function seasonal_UTCI(x::DataFrame)
    p = plot(x.Month, x.UTCI; label="UTCI", color=:red, linewidth=1.5)
    plot!(
        p, x.Month, (x.T2m); label=L"$T_{2m}$", color=:blue, linestyle=:dash, linewidth=1.5
    )
    plot!(p; xlims=(1, 12), xlabel="Month", ylabel="T (°C)", title="Seasonal", grid=true)

    return p
end

function seasonal_Tmrt(x::DataFrame)
    p = plot(x.Month, (x.Tmrt); label=L"$T_{mrt}$", color=:red, linewidth=1.5)
    plot!(
        p, x.Month, (x.T2m); label=L"$T_{2m}$", color=:blue, linestyle=:dash, linewidth=1.5
    )
    plot!(p; xlims=(1, 12), xlabel="Month", ylabel="T (°C)", title="Seasonal", grid=true)

    return p
end

function outdoor_thermal_comfort(TTUrbanDiurnal::DataFrame, TTUrbanSeasonal::DataFrame)
    p1 = diurnal_UTCI(TTUrbanDiurnal)
    p2 = diurnal_Tmrt(TTUrbanDiurnal)
    p3 = seasonal_UTCI(TTUrbanSeasonal)
    p4 = seasonal_Tmrt(TTUrbanSeasonal)

    fig = plot(p1, p2, p3, p4; layout=(2, 2), size=(800, 600))
    return fig
end

function surface_temperature_diurnal(
    x::DataFrame, CFrimp::FT, CFrveg::FT, CFgimp::FT, CFgbare::FT, CFgveg::FT, CFtree::FT
) where {FT<:AbstractFloat}
    p = plot(
        x.Hour, (x.TRoofImp) .* CFrimp; label=L"$T_{roof,imp}$", color=:red, linewidth=1.5
    )
    plot!(
        p,
        x.Hour,
        (x.TRoofVeg) .* CFrveg;
        label=L"$T_{roof,veg}$",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p,
        x.Hour,
        (x.TGroundImp) .* CFgimp;
        label=L"$T_{ground,imp}$",
        color=:red,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(
        p,
        x.Hour,
        (x.TGroundBare) .* CFgbare;
        label=L"$T_{ground,bare}$",
        color=:yellow,
        linewidth=1.5,
    )
    plot!(
        p,
        x.Hour,
        (x.TGroundVeg) .* CFgveg;
        label=L"$T_{ground,veg}$",
        color=:blue,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(p, x.Hour, (x.TTree) .* CFtree; label=L"$T_{tree}$", color=:green, linewidth=1.5)
    plot!(p, x.Hour, (x.TWallSun); label=L"$T_{wall,sun}$", color=:magenta, linewidth=1.5)
    plot!(p, x.Hour, (x.TWallShade); label=L"$T_{wall,shade}$", color=:cyan, linewidth=1.5)
    plot!(
        p;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="T (°C)",
        title="Surface temperature\nDiurnal",
        legend=:bottom,
        legend_columns=2,
        grid=true,
    )

    return p
end

function air_temperature_diurnal(x::DataFrame)
    p = plot(x.Hour, (x.T2m); label=L"$T_{air,2m}$", color=:blue, linewidth=1.5)
    plot!(
        p,
        x.Hour,
        (x.TCanyon);
        label=L"$T_{air,can}$",
        color=:blue,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(p, x.Hour, (x.Tbin); label=L"$T_{b,indoors}$", color=:cyan, linewidth=1.5)
    plot!(p, x.Hour, (x.Tatm); label=L"$T_{air,atm}$", color=:black, linewidth=1.5)
    plot!(
        p;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="T (°C)",
        title="Air temperature\nDiurnal",
        legend=:bottom,
        legend_columns=2,
        grid=true,
    )

    return p
end

function relative_humidity_diurnal(x::DataFrame)
    p = plot(x.Hour, x.RH2m; label=L"$RH_{2m}$", color=:green, linewidth=1.5)
    plot!(p, x.Hour, x.RHbin; label=L"$RH_{b,indoors}$", color=:yellow, linewidth=1.5)
    plot!(p, x.Hour, x.RHatm; label=L"$RH_{atm}$", color=:black, linewidth=1.5)
    plot!(
        p;
        xlims=(0, 23),
        xlabel="hour",
        ylabel="Relative humidity (%)",
        title="Relative humidity\nDiurnal",
        legend=:bottom,
        legend_columns=2,
        grid=true,
    )

    return p
end

function surface_temperature_seasonal(
    x::DataFrame, CFrimp::FT, CFrveg::FT, CFgimp::FT, CFgbare::FT, CFgveg::FT, CFtree::FT
) where {FT<:AbstractFloat}
    p = plot(
        x.Month, (x.TRoofImp) .* CFrimp; label=L"$T_{roof,imp}$", color=:red, linewidth=1.5
    )
    plot!(
        p,
        x.Month,
        (x.TRoofVeg) .* CFrveg;
        label=L"$T_{roof,veg}$",
        color=:blue,
        linewidth=1.5,
    )
    plot!(
        p,
        x.Month,
        (x.TGroundImp) .* CFgimp;
        label=L"$T_{ground,imp}$",
        color=:red,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(
        p,
        x.Month,
        (x.TGroundBare) .* CFgbare;
        label=L"$T_{ground,bare}$",
        color=:yellow,
        linewidth=1.5,
    )
    plot!(
        p,
        x.Month,
        (x.TGroundVeg) .* CFgveg;
        label=L"$T_{ground,veg}$",
        color=:blue,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(p, x.Month, (x.TTree) .* CFtree; label=L"$T_{tree}$", color=:green, linewidth=1.5)
    plot!(p, x.Month, (x.TWallSun); label=L"$T_{wall,sun}$", color=:magenta, linewidth=1.5)
    plot!(p, x.Month, (x.TWallShade); label=L"$T_{wall,shade}$", color=:cyan, linewidth=1.5)
    plot!(p; xlims=(1, 12), xlabel="Month", ylabel="T (°C)", title="Seasonal", grid=true)

    return p
end

function air_temperature_seasonal(x::DataFrame)
    p = plot(x.Month, (x.T2m); label=L"$T_{air,2m}$", color=:blue, linewidth=1.5)
    plot!(
        p,
        x.Month,
        (x.TCanyon);
        label=L"$T_{air,can}$",
        color=:blue,
        linestyle=:dash,
        linewidth=1.5,
    )
    plot!(p, x.Month, (x.Tbin); label=L"$T_{b,indoors}$", color=:cyan, linewidth=1.5)
    plot!(p, x.Month, (x.Tatm); label=L"$T_{air,atm}$", color=:black, linewidth=1.5)
    plot!(p; xlims=(1, 12), xlabel="Month", ylabel="T (°C)", title="Seasonal", grid=true)

    return p
end

function relative_humidity_seasonal(x::DataFrame)
    p = plot(x.Month, x.RH2m; label=L"$RH_{2m}$", color=:green, linewidth=1.5)
    plot!(p, x.Month, x.RHbin; label=L"$RH_{b,indoors}$", color=:yellow, linewidth=1.5)
    plot!(p, x.Month, x.RHatm; label=L"$RH_{atm}$", color=:black, linewidth=1.5)
    plot!(
        p;
        xlims=(1, 12),
        xlabel="Month",
        ylabel="Relative humidity (%)",
        title="Seasonal",
        grid=true,
    )

    return p
end

function surface_temperatures(
    TTUrbanDiurnal::DataFrame, TTUrbanSeasonal::DataFrame, model::Model{FT}
) where {FT<:AbstractFloat}
    CFrimp = model.parameters.surfacefractions.roof.fimp > 0 ? one(FT) : FT(NaN)
    CFrveg = model.parameters.surfacefractions.roof.fveg > 0 ? one(FT) : FT(NaN)
    CFgimp = model.parameters.surfacefractions.ground.fimp > 0 ? one(FT) : FT(NaN)
    CFgbare = model.parameters.surfacefractions.ground.fbare > 0 ? one(FT) : FT(NaN)
    CFgveg = model.parameters.surfacefractions.ground.fveg > 0 ? one(FT) : FT(NaN)
    CFtree = model.parameters.urbangeometry.trees ? one(FT) : FT(NaN)

    p1 = surface_temperature_diurnal(
        TTUrbanDiurnal, CFrimp, CFrveg, CFgimp, CFgbare, CFgveg, CFtree
    )
    p2 = air_temperature_diurnal(TTUrbanDiurnal)
    p3 = relative_humidity_diurnal(TTUrbanDiurnal)
    p4 = surface_temperature_seasonal(
        TTUrbanSeasonal, CFrimp, CFrveg, CFgimp, CFgbare, CFgveg, CFtree
    )
    p5 = air_temperature_seasonal(TTUrbanSeasonal)
    p6 = relative_humidity_seasonal(TTUrbanSeasonal)
    fig = plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(1500, 900))

    return fig
end
