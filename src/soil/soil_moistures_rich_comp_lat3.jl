"""
    soil_moistures_rich_comp_lat3(
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
        Cimp::FT,
        Cbare::FT,
        Cveg::FT,
        fimp::FT,
        fbare::FT,
        fveg::FT,
        Wcan::FT
    ) where {FT<:AbstractFloat}

Calculate lateral soil moisture redistribution between three soil columns (impervious, bare, vegetated).

# Returns
- `dVlat`: Change rates of soil water volume per unit area for the three columns [mm/h]
"""
function soil_moistures_rich_comp_lat3(
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
    Cimp::FT,
    Cbare::FT,
    Cveg::FT,
    fimp::FT,
    fbare::FT,
    fveg::FT,
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
    a = FT(15)        # Same conductivity in horizontal and vertical
    dxsoil = FT(1000) # [mm] = 1 [m]

    # Hydraulic conductivities and soil water potentials
    Kf_gimp = Ko[1]  # [mm/h]
    Kf_gbare = Ko[2] # [mm/h]
    Kf_gveg = Ko[3]  # [mm/h]

    Psi_soil_gimp = Po[1]  # [mm]
    Psi_soil_gbare = Po[2] # [mm]
    Psi_soil_gveg = Po[3]  # [mm]

    # Calculate lateral fluxes between columns
    Qlat_bare2imp = -a * (Kf_gbare + Kf_gimp)/2 * (Psi_soil_gbare - Psi_soil_gimp)/dxsoil
    Qlat_veg2imp = -a * (Kf_gveg + Kf_gimp)/2 * (Psi_soil_gveg - Psi_soil_gimp)/dxsoil
    Qlat_veg2bare = -a * (Kf_gbare + Kf_gveg)/2 * (Psi_soil_gveg - Psi_soil_gbare)/dxsoil
    Qlat_imp2bare = -a * (Kf_gbare + Kf_gimp)/2 * (Psi_soil_gimp - Psi_soil_gbare)/dxsoil
    Qlat_bare2veg = -a * (Kf_gbare + Kf_gveg)/2 * (Psi_soil_gbare - Psi_soil_gveg)/dxsoil
    Qlat_imp2veg = -a * (Kf_gveg + Kf_gimp)/2 * (Psi_soil_gimp - Psi_soil_gveg)/dxsoil

    # Calculate transmissivity
    Tveg2imp = Qlat_veg2imp * dz
    Tbare2imp = Qlat_bare2imp * dz
    Tveg2bare = Qlat_veg2bare * dz
    Timp2bare = Qlat_imp2bare * dz
    Tbare2veg = Qlat_bare2veg * dz
    Timp2veg = Qlat_imp2veg * dz

    # Verify conservation
    T_totflux = sum([Tbare2veg, Tveg2bare, Tbare2imp, Timp2bare, Tveg2imp, Timp2veg])
    if !isapprox(T_totflux, 0; atol=1e-10)
        @warn "The lateral transmissivities do not add up to 0"
    end

    # Calculate incoming fluxes
    Qin_veg2imp = fimp > 0 ? Tveg2imp/(fimp * 1000 * Wcan) * Cveg * Cimp : 0.0
    Qin_bare2imp = fimp > 0 ? Tbare2imp/(fimp * 1000 * Wcan) * Cbare * Cimp : 0.0

    Qin_veg2bare = fbare > 0 ? Tveg2bare/(fbare * 1000 * Wcan) * Cveg * Cbare : 0.0
    Qin_imp2bare = fbare > 0 ? Timp2bare/(fbare * 1000 * Wcan) * Cimp * Cbare : 0.0

    Qin_bare2veg = fveg > 0 ? Tbare2veg/(fveg * 1000 * Wcan) * Cbare * Cveg : 0.0
    Qin_imp2veg = fveg > 0 ? Timp2veg/(fveg * 1000 * Wcan) * Cimp * Cveg : 0.0

    # Sum fluxes for each column
    dVlat = [
        Qin_veg2imp + Qin_bare2imp,    # impervious
        Qin_veg2bare + Qin_imp2bare,   # bare
        Qin_bare2veg + Qin_imp2veg,     # vegetated
    ]

    if any(isnan.(dVlat))
        @warn "NaN values in dVlat"
    end

    return dVlat
end
