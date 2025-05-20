"""
    water_soil(
        Otm1::Vector{FT},
        f::FT,
        TE_H::Vector{FT},
        TE_L::Vector{FT},
        E_soil::FT,
        Qlat_in::Vector{FT},
        dth::FT,
        Pcla::FT,
        Psan::FT,
        Porg::FT,
        Kfc::FT,
        Phy::FT,
        SPAR::Int,
        Kbot::FT,
        CASE_ROOT_H::Int,
        CASE_ROOT_L::Int,
        ZR95_H::Vector{FT},
        ZR95_L::Vector{FT},
        ZR50_H::Vector{FT},
        ZR50_L::Vector{FT},
        ZRmax_H::Vector{FT},
        ZRmax_L::Vector{FT},
        Rrootl_H::Vector{FT},
        Rrootl_L::Vector{FT},
        PsiL50_H::Vector{FT},
        PsiL50_L::Vector{FT},
        PsiX50_H::Vector{FT},
        PsiX50_L::Vector{FT},
        Zs::Vector{FT},
        row::FT
    ) where {FT<:AbstractFloat}

Calculate soil water dynamics and hydrological fluxes.

# Arguments
- `Otm1`: Previous timestep soil moisture [-]
- `f`: Infiltration rate [mm/h]
- `TE_H`: High vegetation transpiration [kg/m²s]
- `TE_L`: Low vegetation transpiration [kg/m²s]
- `E_soil`: Soil evaporation [kg/m²s]
- `Qlat_in`: Lateral water flux [mm/dth]
- `dth`: Time step [h]
- `Pcla`: Clay fraction [-]
- `Psan`: Sand fraction [-]
- `Porg`: Organic matter fraction [-]
- `Kfc`: Hydraulic conductivity at field capacity [mm/h]
- `Phy`: Soil water potential at hygroscopic point [kPa]
- `SPAR`: Soil parameterization choice [-]
- `Kbot`: Bottom boundary conductivity [mm/h]
- `CASE_ROOT_H`: High vegetation root distribution type [-]
- `CASE_ROOT_L`: Low vegetation root distribution type [-]
- `ZR95_H`: 95th percentile root depth for high vegetation [mm]
- `ZR95_L`: 95th percentile root depth for low vegetation [mm]
- `ZR50_H`: 50th percentile root depth for high vegetation [mm]
- `ZR50_L`: 50th percentile root depth for low vegetation [mm]
- `ZRmax_H`: Maximum root depth for high vegetation [mm]
- `ZRmax_L`: Maximum root depth for low vegetation [mm]
- `Rrootl_H`: Root length density for high vegetation [mm/mm³]
- `Rrootl_L`: Root length density for low vegetation [mm/mm³]
- `PsiL50_H`: Leaf water potential at 50% loss for high vegetation [MPa]
- `PsiL50_L`: Leaf water potential at 50% loss for low vegetation [MPa]
- `PsiX50_H`: Xylem water potential at 50% loss for high vegetation [MPa]
- `PsiX50_L`: Xylem water potential at 50% loss for low vegetation [MPa]
- `Zs`: Soil layer depths [mm]
- `row`: Water density [kg/m³]
"""
function water_soil(
    Otm1::Vector{FT},
    f::FT,
    TE_H::FT,
    TE_L::FT,
    E_soil::FT,
    Qlat_in::Vector{FT},
    dth::FT,
    Pcla::FT,
    Psan::FT,
    Porg::FT,
    Kfc::FT,
    Phy::FT,
    SPAR::Int,
    Kbot::FT,
    CASE_ROOT_H::Int,
    CASE_ROOT_L::Int,
    ZR95_H::Vector{FT},
    ZR95_L::Vector{FT},
    ZR50_H::Vector{FT},
    ZR50_L::Vector{FT},
    ZRmax_H::Vector{FT},
    ZRmax_L::Vector{FT},
    Rrootl_H::Vector{FT},
    Rrootl_L::Vector{FT},
    PsiL50_H::Vector{FT},
    PsiL50_L::Vector{FT},
    PsiX50_H::Vector{FT},
    PsiX50_L::Vector{FT},
    Zs::Vector{FT},
    row::FT;
    atol=FT(0.05),
) where {FT<:AbstractFloat}
    # Get soil parameters
    Zs, dz, ms, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, EvL_Zs, Inf_Zs, RfH_Zs, RfL_Zs, _, Kbot, Slo_pot, Dz, aR, aTop = Soil.soil_parameters_total(
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

    # Convert units
    E_soil = E_soil * 3600 * 1000 / row  # [mm/h]
    TE_L = TE_L * 3600 * 1000 / row      # [mm/h]
    TE_H = TE_H * 3600 * 1000 / row      # [mm/h]

    # Distributed sinks
    E_soil_dis = E_soil .* EvL_Zs    # [mm/h] Evaporation from bare soil
    TE_dis_H = TE_H .* vec(RfH_Zs)        # [mm/h] Root water uptake
    TE_dis_L = TE_L .* vec(RfL_Zs)        # [mm/h] Root water uptake

    # Bottom leakage [mm/h]
    Lk = Soil.leakage_bottom(Otm1, Ks_Zs, Osat, Ohy, L, nVG, Kbot, ms, SPAR)

    # Convert lateral inflow to hourly rate
    Qlat_in = Qlat_in / dth  # [mm/h]

    # Initial condition - water volume excluding residual content
    V0 = (Otm1 .- Ohy) .* dz  # [mm]

    # Solve Richards equation
    T_SPAN = (zero(FT), dth)
    ISeep = ones(FT, ms)

    # ODE settings
    maxstep = dth

    # Solve ODE system
    prob = ODEProblem(
        (V, p, t) -> Soil.soil_moistures_rich_comp(
            V,
            Lk,
            f,
            E_soil_dis,
            TE_dis_H,
            TE_dis_L,
            Qlat_in,
            ISeep,
            SPAR,
            Osat,
            Ohy,
            O33,
            dz,
            Ks_Zs,
            Dz,
            ms,
            L,
            Pe,
            aR,
            aTop,
            alpVG,
            nVG,
            one(FT),
            zero(FT),
            zero(FT),
        ),
        V0,
        T_SPAN,
    )

    sol = solve(prob, Rosenbrock23(; autodiff=AutoFiniteDiff()); abstol=atol, dtmax=maxstep)
    V = sol.u[end]

    if any(isnan, V)
        @warn "NaN values in the Volumes"
        return nothing
    end

    # Update soil moisture and calculate outputs
    O, _, _, OS, Psi_s_H, Psi_s_L, _, _, Exwat_H, Exwat_L, Rd, WTR, _, _, _ = Soil.soil_water_multilayer(
        V,
        Zs,
        dz,
        ms,
        Osat,
        Ohy,
        nVG,
        alpVG,
        Ks_Zs,
        L,
        Pe,
        O33,
        SPAR,
        EvL_Zs,
        Inf_Zs,
        RfH_Zs,
        RfL_Zs,
        Rrootl_H,
        Rrootl_L,
        PsiL50_H,
        PsiL50_L,
        PsiX50_H,
        PsiX50_L,
    )

    # Volume correction for water table rise and runoff
    V[1] = V[1] + WTR[2] - Rd
    V[2:(end - 1)] = V[2:(end - 1)] .+ (WTR[3:end] .- WTR[2:(end - 1)])
    V[end] = V[end] - WTR[end]

    # Volume correction for negative values
    if any(<(0), V)
        V, TE_dis_H, TE_dis_L, E_soil, Lk = Soil.volume_correction(
            V, EvL_Zs, RfH_Zs, RfL_Zs, dth*E_soil, dth*TE_dis_H, dth*TE_dis_L, dth*Lk
        )
        TE_H = sum(TE_dis_H) / dth  # [mm/h]
        TE_L = sum(TE_dis_L) / dth  # [mm/h]
        E_soil = E_soil / dth       # [mm/h]
        Lk = Lk / dth              # [mm/h]
    end

    # Calculate potentials and conductivities
    Ko, Psi_soil = Soil.conductivity_suction(
        SPAR, Ks_Zs, Osat, Ohy, L, Pe, O33, alpVG, nVG, O
    )

    # Volume balance check
    dV_dt = sum(V) - sum(V0)
    WBalance_soil =
        dth*f + dth*sum(Qlat_in) - dth*Lk - dth*sum(TE_H) - dth*sum(TE_L) - dth*E_soil -
        Rd - dV_dt

    # Convert units back
    E_soil = E_soil / (3600 * 1000 / row)  # [kg/m²s]
    TE_L = TE_L / (3600 * 1000 / row)      # [kg/m²s]
    TE_H = TE_H / (3600 * 1000 / row)      # [kg/m²s]

    return V,
    O,
    OS,
    Lk,
    Psi_s_H,
    Psi_s_L,
    vec(Exwat_H),
    vec(Exwat_L),
    Rd,
    TE_L,
    TE_H,
    E_soil,
    dV_dt,
    WBalance_soil,
    Psi_soil,
    Ko
end
