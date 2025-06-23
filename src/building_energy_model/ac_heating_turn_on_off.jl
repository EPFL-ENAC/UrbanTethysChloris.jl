"""
    ac_heating_turn_on_off(
        ParHVAC::NamedTuple,
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
- `ParHVAC`: Updated HVAC parameters
- `ParHVACorig`: Original HVAC parameters
"""
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
