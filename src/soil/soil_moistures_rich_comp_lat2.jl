"""
    soil_moistures_rich_comp_lat2(
        Vlat::Vector{FT},
        dz::FT,
        SPAR::Int,
        Ks::FT,
        Osat::FT,
        Ohy::FT,
        L::FT,
        Pe::FT,
        O33::FT,
        alpVG::FT,
        nVG::FT,
        C1::FT,
        C2::FT,
        f1::FT,
        f2::FT,
        Wcan::FT
    ) where {FT<:AbstractFloat}

Calculate lateral soil moisture redistribution between soil columns.

# Arguments
- `Vlat`: Vector of soil water volumes per unit area for both columns [mm]
- `dz`: Layer thickness [mm]
- `SPAR`: Soil parameterization choice (1: Van-Genuchten, 2: Saxton-Rawls)
- `Ks`: Saturated hydraulic conductivity [mm/h]
- `Osat`: Saturated water content [m³/m³]
- `Ohy`: Hygroscopic water content [m³/m³]
- `L`: Lambda parameter (for Saxton-Rawls model) [-]
- `Pe`: Air entry pressure [kPa]
- `O33`: Water content at -33 kPa [m³/m³]
- `alpVG`: Van Genuchten α parameter [1/mm]
- `nVG`: Van Genuchten n parameter [-]
- `C1`: Column 1 coefficient [-]
- `C2`: Column 2 coefficient [-]
- `f1`: Fraction for column 1 [-]
- `f2`: Fraction for column 2 [-]
- `Wcan`: Canopy width [m]

# Returns
- `dVlat`: Change rates of soil water volume per unit area for both columns [mm/h]

# Notes
- If `f1` or `f2` is zero, the corresponding flux is set to zero.
"""
function soil_moistures_rich_comp_lat2(
    Vlat::Vector{FT},
    dz::FT,
    SPAR::Int,
    Ks::FT,
    Osat::FT,
    Ohy::FT,
    L::FT,
    Pe::FT,
    O33::FT,
    alpVG::FT,
    nVG::FT,
    C1::FT,
    C2::FT,
    f1::FT,
    f2::FT,
    Wcan::FT,
) where {FT<:AbstractFloat}
    # Calculate water content from volume
    Olat = Vlat ./ dz .+ Ohy

    # Bound water content between hygroscopic and near-saturation
    I1 = Olat .>= Osat - 1e-5
    I2 = Olat .<= Ohy + 1e-5
    Olat[I1] .= Osat - 1e-5
    Olat[I2] .= Ohy + 1e-5

    # Get conductivity and suction using existing function
    results = conductivity_suction.(SPAR, Ks, Osat, Ohy, L, Pe, O33, alpVG, nVG, Olat)
    Ko = [x[1] for x in results]
    Po = [x[2] for x in results]

    # Lateral water redistribution parameters
    a = 15.0        # Same conductivity in horizontal and vertical
    dxsoil = 1000.0 # [mm] = 1 [m]

    # Calculate lateral fluxes between columns
    Qlat_1to2 = -a * (Ko[1] + Ko[2])/2 * (Po[1] - Po[2])/dxsoil
    Qlat_2to1 = -a * (Ko[1] + Ko[2])/2 * (Po[2] - Po[1])/dxsoil

    # Calculate transmissivity
    T_1to2 = Qlat_1to2 * dz
    T_2to1 = Qlat_2to1 * dz

    # Verify conservation
    if T_1to2 + T_2to1 ≉ 0.0
        @warn "The lateral transmissivities do not add up to 0"
    end

    # Calculate incoming fluxes
    Qin_1to2 = f2 > 0 ? T_1to2/(f2 * 1000 * Wcan) * C1 * C2 : 0.0
    Qin_2to1 = f1 > 0 ? T_2to1/(f1 * 1000 * Wcan) * C2 * C1 : 0.0

    return [Qin_2to1, Qin_1to2]
end
