"""
    ac_heating_module(
        AC_on::Bool,
        Heat_on::Bool,
        AC_onCool::Bool,
        AC_onDehum::Bool,
        ParHVAC::NamedTuple,
        HbuildIn::FT,
        Hvent::FT,
        Hequip::FT,
        Hpeople::FT,
        dSH_air::FT,
        LEvent::FT,
        LEequip::FT,
        LEpeople::FT,
        dSLE_air::FT
    ) where {FT<:AbstractFloat}

AC & heating module that removes or adds sensible and latent heat as needed to maintain air temperature
and humidity. For AC, both temperature and humidity are controlled. For heating, only temperature is controlled.

# Arguments
- `AC_on`: Master AC switch (on/off)
- `Heat_on`: Heating switch (on/off)
- `AC_onCool`: AC cooling switch (on/off)
- `AC_onDehum`: AC dehumidification switch (on/off)
- `ParHVAC`: HVAC parameters
- `HbuildIn`: Building internal sensible heat [W]
- `Hvent`: Ventilation sensible heat [W]
- `Hequip`: Equipment sensible heat [W]
- `Hpeople`: People sensible heat [W]
- `dSH_air`: Change in sensible heat storage [W]
- `LEvent`: Ventilation latent heat [W]
- `LEequip`: Equipment latent heat [W]
- `LEpeople`: People latent heat [W]
- `dSLE_air`: Change in latent heat storage [W]

# Returns
- `AC_on::Bool`: Master AC switch status
- `AC_onCool::Bool`: AC cooling switch status
- `AC_onDehum::Bool`: AC dehumidification switch status
- `Heat_on::Bool`: Heating switch status
- `H_AC_Heat::FT`: Sensible heat removed/added by HVAC [W]
- `LE_AC_Heat::FT`: Latent heat removed/added by HVAC [W]
"""
function ac_heating_module(
    AC_on::Bool,
    Heat_on::Bool,
    AC_onCool::Bool,
    AC_onDehum::Bool,
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    HbuildIn::FT,
    Hvent::FT,
    Hequip::FT,
    Hpeople::FT,
    dSH_air::FT,
    LEvent::FT,
    LEequip::FT,
    LEpeople::FT,
    dSLE_air::FT,
) where {FT<:AbstractFloat}

    # Initialize outputs
    H_AC_Heat = zero(FT)
    LE_AC_Heat = zero(FT)

    # Cooling mode
    if AC_on
        if AC_onCool && AC_onDehum
            H_AC_Heat = HbuildIn + Hvent + Hequip + Hpeople - dSH_air
            LE_AC_Heat = LEvent + LEequip + LEpeople - dSLE_air

        elseif AC_onCool && !AC_onDehum
            H_AC_Heat = HbuildIn + Hvent + Hequip + Hpeople - dSH_air
            LE_AC_Heat = zero(FT)

        elseif !AC_onCool && AC_onDehum
            H_AC_Heat = zero(FT)
            LE_AC_Heat = LEvent + LEequip + LEpeople - dSLE_air
        end

        # Heating mode
    elseif Heat_on
        H_AC_Heat = HbuildIn + Hvent + Hequip + Hpeople - dSH_air
        LE_AC_Heat = zero(FT)
    end

    return AC_on::Bool,
    AC_onCool::Bool, AC_onDehum::Bool, Heat_on::Bool, H_AC_Heat::FT,
    LE_AC_Heat::FT
end
