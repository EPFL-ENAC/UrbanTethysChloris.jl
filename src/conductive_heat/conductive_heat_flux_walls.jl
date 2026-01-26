"""
    conductive_heat_flux_walls(
        TemperatureC::Vector{FT},
        TemperatureB::Vector{FT},
        TempVec_ittm::NamedTuple,
        TempVecB_ittm::NamedTuple,
        Anthropogenic::NamedTuple,
        ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
        WallLayers::NamedTuple,
        ParCalculation::NamedTuple,
        type::Bool,
        ParWindows::ModelComponents.Parameters.WindowParameters{FT},
        BEM_on::Bool
    ) where {FT<:AbstractFloat}

Calculate conductive heat flux through walls.

# Arguments
- `TemperatureC`: Canyon temperatures vector
- `TemperatureB`: Building temperatures vector
- `TempVec_ittm`: Temperature vectors at previous time step
- `TempVecB_ittm`: Building temperature vectors at previous time step
- `Anthropogenic`: Anthropogenic parameters
- `ParThermalWall`: Thermal parameters for wall
- `WallLayers`: Wall layer parameters
- `ParCalculation`: Calculation parameters
- `type`: 1 for sunlit wall, 0 for shaded wall
- `ParWindows`: Window parameters
- `BEM_on`: Building Energy Model switch

# Returns
- `G1::FT`: Heat flux from surface to concrete interior [W/m²]
- `G2::FT`: Heat flux from concrete interior to building interior [W/m²]
- `dS::FT`: Energy storage in the wall [J/m²]
"""
function conductive_heat_flux_walls_2(
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    TempVecB_ittm::ModelComponents.ModelVariables.TempVecB{FT},
    Anthropogenic::ModelComponents.ForcingInputs.AnthropogenicInputs{FT,0},
    ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    WallLayers::NamedTuple,
    ParCalculation::NamedTuple,
    type::Bool,
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Extract parameters based on wall type
    if type  # sunlit wall
        Ts = TemperatureC[4]
        Tb = BEM_on ? TemperatureB[2] : Anthropogenic.Tb
        Tint = TemperatureC[7]
        Ts_tm1 = TempVec_ittm.TWallSun
        Tb_tm1 = TempVecB_ittm.Tinwallsun
        Tint_tm1 = TempVec_ittm.TWallIntSun
    elseif !type  # shaded wall
        Ts = TemperatureC[5]
        Tb = BEM_on ? TemperatureB[3] : Anthropogenic.Tb
        Tint = TemperatureC[8]
        Ts_tm1 = TempVec_ittm.TWallShade
        Tb_tm1 = TempVecB_ittm.Tinwallshd
        Tint_tm1 = TempVec_ittm.TWallIntShade
    end

    # Extract common parameters
    lan_dry1 = ParThermalWall.lan_dry
    lan_dry2 = ParThermalWall.lan_dry
    dz1 = WallLayers.dz1_wall
    dz2 = WallLayers.dz2_wall
    cv_s1 = ParThermalWall.cv_s
    cv_s2 = ParThermalWall.cv_s
    dts = ParCalculation.dts

    # Window parameters
    UvalueWindow = ParWindows.Uvalue
    cv_glass = ParWindows.cv_glass
    dzglass = ParWindows.dztot
    RatioWindowToWall = min(ParWindows.GlazingRatio, FT(0.95))

    # Compute wall heat fluxes
    G1wall = lan_dry1 * (Ts - Tint) / dz1
    G2wall = lan_dry2 * (Tint - Tb) / dz2
    dSwall = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

    # Compute window heat fluxes
    G1window = UvalueWindow * (Ts - Tb)
    G2window = UvalueWindow * (Ts - Tb)
    dSwindow = cv_glass * dzglass / dts * ((Ts + Tb) / 2 - (Ts_tm1 + Tb_tm1) / 2)

    # Average fluxes based on window ratio when BEM is active
    if BEM_on
        G1 = (1 - RatioWindowToWall) * G1wall + RatioWindowToWall * G1window
        G2 = (1 - RatioWindowToWall) * G2wall + RatioWindowToWall * G2window
        dS = (1 - RatioWindowToWall) * dSwall + RatioWindowToWall * dSwindow
    else
        G1 = G1wall
        G2 = G2wall
        dS = dSwall
    end

    # @infiltrate
    return G1, G2, dS
end
function conductive_heat_flux_walls(
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    TempVecB_ittm::ModelComponents.ModelVariables.TempVecB{FT},
    Anthropogenic::ModelComponents.ForcingInputs.AnthropogenicInputs{FT,0},
    ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    WallLayers::NamedTuple,
    ParCalculation::NamedTuple,
    type::Bool,
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Extract parameters based on wall type
    if type == 1  # sunlit wall
        Ts = TemperatureC[4]
        Tb = BEM_on ? TemperatureB[2] : Anthropogenic.Tb
        Tint = TemperatureC[7]
        Ts_tm1 = TempVec_ittm.TWallSun
        Tb_tm1 = TempVecB_ittm.Tinwallsun
        Tint_tm1 = TempVec_ittm.TWallIntSun
    elseif type == 0  # shaded wall
        Ts = TemperatureC[5]
        Tb = BEM_on ? TemperatureB[3] : Anthropogenic.Tb
        Tint = TemperatureC[8]
        Ts_tm1 = TempVec_ittm.TWallShade
        Tb_tm1 = TempVecB_ittm.Tinwallshd
        Tint_tm1 = TempVec_ittm.TWallIntShade
    else
        error("type must be 0 (shaded) or 1 (sunlit)")
    end

    # Extract common parameters
    lan_dry1 = ParThermalWall.lan_dry
    lan_dry2 = ParThermalWall.lan_dry
    dz1 = WallLayers.dz1_wall
    dz2 = WallLayers.dz2_wall
    cv_s1 = ParThermalWall.cv_s
    cv_s2 = ParThermalWall.cv_s
    dts = ParCalculation.dts

    # Window parameters
    UvalueWindow = ParWindows.Uvalue
    cv_glass = ParWindows.cv_glass
    dzglass = ParWindows.dztot
    RatioWindowToWall = min(ParWindows.GlazingRatio, FT(0.95))

    # Compute wall heat fluxes
    G1wall = lan_dry1 * (Ts - Tint) / dz1
    G2wall = lan_dry2 * (Tint - Tb) / dz2
    dSwall = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

    # Compute window heat fluxes
    G1window = UvalueWindow * (Ts - Tb)
    G2window = UvalueWindow * (Ts - Tb)
    dSwindow = cv_glass * dzglass / dts * ((Ts + Tb) / 2 - (Ts_tm1 + Tb_tm1) / 2)

    # Average fluxes based on window ratio when BEM is active
    if BEM_on
        G1 = (1 - RatioWindowToWall) * G1wall + RatioWindowToWall * G1window
        G2 = (1 - RatioWindowToWall) * G2wall + RatioWindowToWall * G2window
        dS = (1 - RatioWindowToWall) * dSwall + RatioWindowToWall * dSwindow
    else
        G1 = G1wall
        G2 = G2wall
        dS = dSwall
    end

    return G1, G2, dS
end
function conductive_heat_flux_walls(
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVec_ittm::NamedTuple,
    TempVecB_ittm::NamedTuple,
    Anthropogenic::NamedTuple,
    ParThermalWall::ModelComponents.Parameters.LocationSpecificThermalProperties{FT},
    WallLayers::NamedTuple,
    ParCalculation::NamedTuple,
    type::Bool,
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Extract parameters based on wall type
    if type == 1  # sunlit wall
        Ts = TemperatureC[4]
        Tb = BEM_on ? TemperatureB[2] : Anthropogenic.Tb
        Tint = TemperatureC[7]
        Ts_tm1 = TempVec_ittm.TWallSun
        Tb_tm1 = TempVecB_ittm.Tinwallsun
        Tint_tm1 = TempVec_ittm.TWallIntSun
    elseif type == 0  # shaded wall
        Ts = TemperatureC[5]
        Tb = BEM_on ? TemperatureB[3] : Anthropogenic.Tb
        Tint = TemperatureC[8]
        Ts_tm1 = TempVec_ittm.TWallShade
        Tb_tm1 = TempVecB_ittm.Tinwallshd
        Tint_tm1 = TempVec_ittm.TWallIntShade
    else
        error("type must be 0 (shaded) or 1 (sunlit)")
    end

    # Extract common parameters
    lan_dry1 = ParThermalWall.lan_dry
    lan_dry2 = ParThermalWall.lan_dry
    dz1 = WallLayers.dz1_wall
    dz2 = WallLayers.dz2_wall
    cv_s1 = ParThermalWall.cv_s
    cv_s2 = ParThermalWall.cv_s
    dts = ParCalculation.dts

    # Window parameters
    UvalueWindow = ParWindows.Uvalue
    cv_glass = ParWindows.cv_glass
    dzglass = ParWindows.dztot
    RatioWindowToWall = min(ParWindows.GlazingRatio, FT(0.95))

    # Compute wall heat fluxes
    G1wall = lan_dry1 * (Ts - Tint) / dz1
    G2wall = lan_dry2 * (Tint - Tb) / dz2
    dSwall = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

    # Compute window heat fluxes
    G1window = UvalueWindow * (Ts - Tb)
    G2window = UvalueWindow * (Ts - Tb)
    dSwindow = cv_glass * dzglass / dts * ((Ts + Tb) / 2 - (Ts_tm1 + Tb_tm1) / 2)

    # Average fluxes based on window ratio when BEM is active
    if BEM_on
        G1 = (1 - RatioWindowToWall) * G1wall + RatioWindowToWall * G1window
        G2 = (1 - RatioWindowToWall) * G2wall + RatioWindowToWall * G2window
        dS = (1 - RatioWindowToWall) * dSwall + RatioWindowToWall * dSwindow
    else
        G1 = G1wall
        G2 = G2wall
        dS = dSwall
    end

    return G1, G2, dS
end
