"""
    soil_thermal_properties(Tdp, rsd, lan_dry, lan_s, cv_s, Osat, Ohy, O) where {FT<:AbstractFloat}

Compute soil thermal properties based on soil moisture and composition.

# Arguments
- `Tdp::Vector{FT},`: Dampening temperature [°C]
- `rsd::Vector{FT}`: Density dry soil [kg/m³]
- `lan_dry::Vector{FT}`: Thermal conductivity dry soil [W/m K]
- `lan_s::Vector{FT}`: Thermal conductivity soil solid [W/m K]
- `cv_s::Vector{FT}`: Volumetric heat capacity soil solid [J/m³ K]
- `Osat::Vector{FT}`: Water content at saturation [-]
- `Ohy::Vector{FT}`: Residual/hygroscopic water content [-]
- `O::Vector{FT}`: Soil moisture at previous time step [-]

# Returns
- `lanS::Vector{FT}`: Thermal conductivity of soil [W/m K]
- `cv_Soil::Vector{FT}`: Volumetric heat capacity of soil [J/m³ K]
- `CTt::Vector{FT}`: Total thermal capacity of soil [K m²/J]
"""
function soil_thermal_properties(
    Tdp::Vector{FT},
    rsd::Vector{FT},
    lan_dry::Vector{FT},
    lan_s::Vector{FT},
    cv_s::Vector{FT},
    Osat::Vector{FT},
    Ohy::Vector{FT},
    O::Vector{FT},
) where {FT<:AbstractFloat}
    # Water properties
    row = 1000.0  # [kg/m^3]
    lan_wat = 0.58  # [W/m K] Thermal conductivity water
    lan_ice = 2.29  # [W/m K] Thermal conductivity ice
    cv_w = 4186000.0  # [J/m^3 K] Volumetric heat capacity water

    n = length(O)
    lanS = zeros(n)    # [W/m K] Thermal conductivity Soil
    cv_Soil = zeros(n) # Volumetric heat capacity Soil [J/m^3 K]
    rsoil = zeros(n)   # Soil Density [kg/m^3]
    cs_Soil = zeros(n) # [J/kg K] Specific Heat Soil

    # Handle boundary conditions
    O = clamp.(O, Ohy, Osat)
    Oice = Osat .* (Tdp .< 0)  # Frozen layers

    # Process each soil layer
    for i in 1:n
        if Tdp[i] > 0
            lan_sat = lan_wat^Osat[i] * lan_s[i]^(1-Osat[i])
            Ke = log((O[i] + Oice[i])/Osat[i]) + 1
            Ke = max(Ke, 0)
        else
            Oliq = Osat[i] - Oice[i]
            lan_sat = lan_wat^Osat[i] * lan_s[i]^(1-Osat[i]) * lan_ice^(Osat[i]-Oliq)
            Ke = (O[i] + Oice[i])/Osat[i]
        end

        if O[i]/Osat[i] > 1e-7
            lanS[i] = Ke * lan_sat + (1-Ke) * lan_dry[i]
        else
            lanS[i] = lan_dry[i]
        end

        cv_Soil[i] = cv_s[i] * (1-Osat[i]) + O[i] * cv_w
        rsoil[i] = rsd[i] + (O[i]-Ohy[i]) * row
        cs_Soil[i] = cv_Soil[i] / rsoil[i]
    end

    tau = 86400.0  # [s] time constant
    CTt = 2.0 .* sqrt.(π ./ (lanS .* cs_Soil .* rsoil .* tau))

    return lanS, cv_Soil, CTt
end

function soil_thermal_properties(
    Tdp::FT,
    rsd::Vector{FT},
    lan_dry::Vector{FT},
    lan_s::Vector{FT},
    cv_s::Vector{FT},
    Osat::Vector{FT},
    Ohy::Vector{FT},
    O::Vector{FT},
) where {FT<:AbstractFloat}
    # Water properties
    row = 1000.0  # [kg/m^3]
    lan_wat = 0.58  # [W/m K] Thermal conductivity water
    lan_ice = 2.29  # [W/m K] Thermal conductivity ice
    cv_w = 4186000.0  # [J/m^3 K] Volumetric heat capacity water

    n = length(O)
    lanS = zeros(n)    # [W/m K] Thermal conductivity Soil
    cv_Soil = zeros(n) # Volumetric heat capacity Soil [J/m^3 K]
    rsoil = zeros(n)   # Soil Density [kg/m^3]
    cs_Soil = zeros(n) # [J/kg K] Specific Heat Soil

    # Handle boundary conditions
    O = clamp.(O, Ohy, Osat)
    Oice = Osat .* (Tdp < 0)  # Frozen layers

    # Process each soil layer
    for i in 1:n
        if Tdp > 0
            lan_sat = lan_wat^Osat[i] * lan_s[i]^(1-Osat[i])
            Ke = log((O[i] + Oice[i])/Osat[i]) + 1
            Ke = max(Ke, 0)
        else
            Oliq = Osat[i] - Oice[i]
            lan_sat = lan_wat^Osat[i] * lan_s[i]^(1-Osat[i]) * lan_ice^(Osat[i]-Oliq)
            Ke = (O[i] + Oice[i])/Osat[i]
        end

        if O[i]/Osat[i] > 1e-7
            lanS[i] = Ke * lan_sat + (1-Ke) * lan_dry[i]
        else
            lanS[i] = lan_dry[i]
        end

        cv_Soil[i] = cv_s[i] * (1-Osat[i]) + O[i] * cv_w
        rsoil[i] = rsd[i] + (O[i]-Ohy[i]) * row
        cs_Soil[i] = cv_Soil[i] / rsoil[i]
    end

    tau = 86400.0  # [s] time constant
    CTt = 2.0 .* sqrt.(π ./ (lanS .* cs_Soil .* rsoil .* tau))

    return lanS, cv_Soil, CTt
end
