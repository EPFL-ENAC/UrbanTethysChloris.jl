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

"""
    ac_heating_turn_on_off(
        ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
        TempVecB_ittm::NamedTuple,
        TempVec_ittm::NamedTuple,
        Humidity_ittm::NamedTuple,
        MeteoData::NamedTuple,
        Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        BEM_on::Bool
    ) where {FT<:AbstractFloat}

Turn AC and heating on/off based on temperature and humidity conditions.

# Arguments
- `ParHVAC`: HVAC parameters
- `TempVecB_ittm`: Building temperature vectors at previous time step
- `TempVec_ittm`: Temperature vectors at previous time step
- `Humidity_ittm`: Humidity parameters
- `MeteoData`: Meteorological data
- `Geometry_m`: Urban geometry parameters
- `BEM_on`: Building Energy Model switch

# Returns
- `ParHVAC::ModelComponents.Parameters.HVACParameters{FT}`: Updated HVAC parameters
- `ParHVACorig::ModelComponents.Parameters.HVACParameters{FT}`: Original HVAC parameters
"""
function ac_heating_turn_on_off(model::Model{FT}, BEM_on::Bool) where {FT<:AbstractFloat}
    return ac_heating_turn_on_off(
        model.parameters.building_energy.hvac,
        model.variables.buildingenergymodel.TempVecB,
        model.variables.temperature.tempvec,
        model.variables.humidity.Humidity,
        model.forcing.meteorological,
        model.parameters.urbangeometry,
        BEM_on,
    )
end

function ac_heating_turn_on_off(
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    TempVecB_ittm::ModelComponents.ModelVariables.TempVecB{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    Humidity_ittm::ModelComponents.ModelVariables.Humidity{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Store original parameters
    ParHVACorig = ParHVAC

    # Urban geometry
    Hwall = Geometry_m.Height_canyon    # Building height [m]
    Wroof = Geometry_m.Width_roof       # Roof width [m]

    # Previous timestep conditions
    Tbin_tm1 = TempVecB_ittm.Tbin         # Building internal temperature
    qbin_tm1 = TempVecB_ittm.qbin         # Building internal humidity
    Tintmasstm1 = TempVecB_ittm.Tintmass  # Ground temperature

    # Calculate average surface temperature
    Tsurftm1 =
        (
            Wroof * (TempVecB_ittm.Tinground + TempVecB_ittm.Tceiling) +
            Hwall * (TempVecB_ittm.Tinwallshd + TempVecB_ittm.Tinwallsun)
        ) / (2 * Wroof + 2 * Hwall)

    # Extract HVAC parameters
    ACon = ParHVAC.ACon
    Heatingon = ParHVAC.Heatingon
    TsetpointCooling = ParHVAC.TsetpointCooling
    TsetpointHeating = ParHVAC.TsetpointHeating
    RHsetpointCooling = ParHVAC.RHsetpointCooling / 100
    RHsptHeating = NaN  # Not used in the model

    # Calculate humidity metrics for indoor regulation
    esat_TspCooling =
        611 *
        exp(17.27 * (TsetpointCooling - 273.16) / (237.3 + (TsetpointCooling - 273.16)))
    ea_RHspCooling = RHsetpointCooling * esat_TspCooling
    q_RHspCooling = 0.622 * ea_RHspCooling / (MeteoData.Pre - 0.378 * ea_RHspCooling)

    # Check heating and cooling setpoints
    if TsetpointCooling == TsetpointHeating && ACon && Heatingon
        TsetpointHeating = TsetpointHeating - 0.01
    elseif TsetpointHeating > TsetpointCooling && ACon && Heatingon
        @warn "Cooling set point is lower than heating setpoint, please double check"
    end

    # Initialize cooling and dehumidification parameters
    AC_onCool = ACon ? true : false
    AC_onDehum = ACon ? true : false

    # Automatic switch for AC and heating
    if ACon
        if TempVec_ittm.TCanyon >= TsetpointCooling && Tsurftm1 >= TsetpointCooling
            AC_onCool = true
            AC_onDehum = true
            if Humidity_ittm.CanyonSpecific < q_RHspCooling
                AC_onDehum = false
            end
        else
            ACon = false
            AC_onCool = false
            AC_onDehum = false
        end
    end

    if Heatingon
        Heatingon = TempVec_ittm.TCanyon <= TsetpointHeating && Tsurftm1 <= TsetpointHeating
    end

    # Update HVAC parameters based on BEM status
    if BEM_on
        ParHVAC = ConstructionBase.setproperties(
            ParHVAC,
            (
                ACon=ACon,
                AC_onCool=AC_onCool,
                AC_onDehum=AC_onDehum,
                Heatingon=Heatingon,
                TsetpointCooling=TsetpointCooling,
                TsetpointHeating=TsetpointHeating,
                q_RHspCooling=q_RHspCooling,
                MasterOn=false,
            ),
        )
    else
        ParHVAC = ConstructionBase.setproperties(
            ParHVAC,
            (
                ACon=false,
                AC_onCool=false,
                AC_onDehum=false,
                Heatingon=false,
                TsetpointCooling=TsetpointCooling,
                TsetpointHeating=TsetpointHeating,
                q_RHspCooling=q_RHspCooling,
                MasterOn=false,
            ),
        )
    end

    return ParHVAC, ParHVACorig
end
function ac_heating_turn_on_off(
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    TempVecB_ittm::NamedTuple,
    TempVec_ittm::NamedTuple,
    Humidity_ittm::NamedTuple,
    MeteoData::NamedTuple,
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    BEM_on::Bool,
) where {FT<:AbstractFloat}
    # Store original parameters
    ParHVACorig = ParHVAC

    # Urban geometry
    Hwall = Geometry_m.Height_canyon    # Building height [m]
    Wroof = Geometry_m.Width_roof       # Roof width [m]

    # Previous timestep conditions
    Tbin_tm1 = TempVecB_ittm.Tbin         # Building internal temperature
    qbin_tm1 = TempVecB_ittm.qbin         # Building internal humidity
    Tintmasstm1 = TempVecB_ittm.Tintmass  # Ground temperature

    # Calculate average surface temperature
    Tsurftm1 =
        (
            Wroof * (TempVecB_ittm.Tinground + TempVecB_ittm.Tceiling) +
            Hwall * (TempVecB_ittm.Tinwallshd + TempVecB_ittm.Tinwallsun)
        ) / (2 * Wroof + 2 * Hwall)

    # Extract HVAC parameters
    ACon = ParHVAC.ACon
    Heatingon = ParHVAC.Heatingon
    TsetpointCooling = ParHVAC.TsetpointCooling
    TsetpointHeating = ParHVAC.TsetpointHeating
    RHsetpointCooling = ParHVAC.RHsetpointCooling / 100
    RHsptHeating = NaN  # Not used in the model

    # Calculate humidity metrics for indoor regulation
    esat_TspCooling =
        611 *
        exp(17.27 * (TsetpointCooling - 273.16) / (237.3 + (TsetpointCooling - 273.16)))
    ea_RHspCooling = RHsetpointCooling * esat_TspCooling
    q_RHspCooling = 0.622 * ea_RHspCooling / (MeteoData.Pre - 0.378 * ea_RHspCooling)

    # Check heating and cooling setpoints
    if TsetpointCooling == TsetpointHeating && ACon && Heatingon
        TsetpointHeating = TsetpointHeating - 0.01
    elseif TsetpointHeating > TsetpointCooling && ACon && Heatingon
        @warn "Cooling set point is lower than heating setpoint, please double check"
    end

    # Initialize cooling and dehumidification parameters
    AC_onCool = ACon ? true : false
    AC_onDehum = ACon ? true : false

    # Automatic switch for AC and heating
    if ACon
        if TempVec_ittm.TCanyon >= TsetpointCooling && Tsurftm1 >= TsetpointCooling
            AC_onCool = true
            AC_onDehum = true
            if Humidity_ittm.CanyonSpecific < q_RHspCooling
                AC_onDehum = false
            end
        else
            ACon = false
            AC_onCool = false
            AC_onDehum = false
        end
    end

    if Heatingon
        Heatingon = TempVec_ittm.TCanyon <= TsetpointHeating && Tsurftm1 <= TsetpointHeating
    end

    # Update HVAC parameters based on BEM status
    if BEM_on
        ParHVAC = ConstructionBase.setproperties(
            ParHVAC,
            (
                ACon=ACon,
                AC_onCool=AC_onCool,
                AC_onDehum=AC_onDehum,
                Heatingon=Heatingon,
                TsetpointCooling=TsetpointCooling,
                TsetpointHeating=TsetpointHeating,
                q_RHspCooling=q_RHspCooling,
                MasterOn=false,
            ),
        )
    else
        ParHVAC = ConstructionBase.setproperties(
            ParHVAC,
            (
                ACon=false,
                AC_onCool=false,
                AC_onDehum=false,
                Heatingon=false,
                TsetpointCooling=TsetpointCooling,
                TsetpointHeating=TsetpointHeating,
                q_RHspCooling=q_RHspCooling,
                MasterOn=false,
            ),
        )
    end

    return ParHVAC, ParHVACorig
end
