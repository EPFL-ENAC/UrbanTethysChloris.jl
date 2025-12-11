"""
    eb_solver_building_output(
        TemperatureC::Vector{FT},
        TemperatureB::Vector{FT},
        TempVecB_ittm::NamedTuple,
        TempVec_ittm::NamedTuple,
        Humidity_ittm::NamedTuple,
        MeteoData::NamedTuple,
        SWRinWsun::FT,
        SWRinWshd::FT,
        G2Roof::FT,
        G2WallSun::FT,
        G2WallShade::FT,
        TempDamp_ittm::NamedTuple,
        SWRabs_t::NamedTuple,
        Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
        PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
        ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
        ParCalculation::NamedTuple,
        ParThermalBuildingInt::ModelComponents.Parameters.ThermalBuilding{FT},
        ParWindows::ModelComponents.Parameters.WindowParameters{FT},
        BEM_on::Bool,
        HVACSchedule::NamedTuple,
    ) where {FT<:AbstractFloat}

Simple Building energy model.

# Arguments
- `TemperatureC`: Outdoor temperatures of buildings and air and humidity [K]
- `TemperatureB`: Building internal temperature and air and humidity [K]
- `TempVecB_ittm`: Building internal temperature and humidity from previous time step
- `TempVec_ittm`: Temperature vectors at previous time step
- `Humidity_ittm`: Humidity at previous time step
- `MeteoData`: Atmospheric forcing conditions
- `SWRinWsun`: Incoming shortwave radiation on sunlit wall [W/m²]
- `SWRinWshd`: Incoming shortwave radiation on shaded wall [W/m²]
- `G2Roof`: Conductive heat flux through roof [W/m²]
- `G2WallSun`: Conductive heat flux through sunlit wall [W/m²]
- `G2WallShade`: Conductive heat flux through shaded wall [W/m²]
- `TempDamp_ittm`: Ground dampening temperature from previous time step
- `SWRabs_t`: Absorbed shortwave radiation [W/m²]
- `Geometry_m`: Urban geometry parameters
- `PropOpticalIndoors`: Building internal albedo/emissivities
- `ParHVAC`: HVAC parameters
- `ParCalculation`: Calculation parameters
- `ParThermalBuildingInt`: Thermal parameters for building interior
- `ParWindows`: Window parameters
- `BEM_on`: Building Energy Model switch
- `HVACSchedule`: HVAC operation schedule

# Returns
- `HbuildInt::NamedTuple`: Internal sensible heat fluxes [W/m²]
- `LEbuildInt::NamedTuple`: Internal latent heat fluxes [W/m²]
- `GbuildInt::NamedTuple`: Internal conductive heat fluxes [W/m²]
- `SWRabsB::NamedTuple`: Absorbed shortwave radiation [W/m²]
- `LWRabsB::NamedTuple`: Absorbed longwave radiation [W/m²]
- `Tdpfloor`: Floor dampening temperature [K]
- `WasteHeat`: Waste heat emissions [W/m²]
- `EnergyUse`: Energy consumption [W*h]
- `HumidityBuilding`: Building humidity parameters
- `ParACHeat`: HVAC system operation parameters
- `YBuildInt::Vector{FT}`: Building internal energy balance residuals [W/m²]
"""
function eb_solver_building_output(
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVecB_ittm::ModelComponents.ModelVariables.TempVecB{FT},
    TempVec_ittm::ModelComponents.ModelVariables.TempVec{FT},
    Humidity_ittm::ModelComponents.ModelVariables.Humidity{FT},
    MeteoData::ModelComponents.ForcingInputs.MeteorologicalInputs{FT,0},
    SWRinWsun::FT,
    SWRinWshd::FT,
    G2Roof::FT,
    G2WallSun::FT,
    G2WallShade::FT,
    TempDamp_ittm::ModelComponents.ModelVariables.TempDamp{FT},
    SWRabs_t::Radiation.RadiationFluxes{FT},
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    ParCalculation::NamedTuple,
    ParThermalBuildingInt::ModelComponents.Parameters.ThermalBuilding{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
    HVACSchedule::ModelComponents.ForcingInputs.HVACSchedule{FT,0};
) where {FT<:AbstractFloat}
    # Initialize building interior temperatures and humidity
    Tceiling = TemperatureB[1]      # Temperature ceiling [K]
    Tinwallsun = TemperatureB[2]    # Temperature sunlit wall [K]
    Tinwallshd = TemperatureB[3]    # Temperature shaded wall [K]
    Twindow = TemperatureB[4]       # Temperature windows [K]
    Tinground = TemperatureB[5]     # Temperature ground [K]
    Tintmass = TemperatureB[6]      # Temperature internal mass [K]
    Tbin = TemperatureB[7]          # Air temperature building interior [K]
    qbin = TemperatureB[8]          # Specific humidity air building interior [kg/kg]

    # External canyon air temperature and specific humidity
    Tairout = TemperatureC[9]       # Air temperature at outdoor calculation height [K]
    qairout = TemperatureC[10]      # Specific humidity at outdoor calculation height [kg/kg]

    # Urban geometry [m]
    Hwall = Geometry_m.Height_canyon  # Building height [m]
    Wroof = Geometry_m.Width_roof     # Roof width [m]
    Wcan = Geometry_m.Width_canyon    # Street/canyon width [m]

    # Temperature [K] and specific humidity [kg/kg] conditions from previous time step
    Tbin_tm1 = TempVecB_ittm.Tbin             # Building internal temperature
    qbin_tm1 = TempVecB_ittm.qbin             # Building internal humidity
    Tingroundtm1 = TempVecB_ittm.Tinground    # Ground temperature
    TingroundDamptm1 = TempDamp_ittm.TDampGroundBuild  # Within ground dampening temperature
    Tintmasstm1 = TempVecB_ittm.Tintmass      # Ground temperature

    Tsurftm1 =
        (
            Wroof * (TempVecB_ittm.Tinground + TempVecB_ittm.Tceiling) +
            Hwall * (TempVecB_ittm.Tinwallshd + TempVecB_ittm.Tinwallsun)
        ) / (2 * Wroof + 2 * Hwall)

    # Atmospheric forcing conditions
    Tatm = MeteoData.Tatm  # Air temperature at atmospheric forcing height [K]
    Pre = MeteoData.Pre    # Air pressure at atmospheric forcing height [Pa]
    ea = MeteoData.ea      # Vapor pressure at atmospheric forcing height [Pa]

    # Albedo [-]
    abc = PropOpticalIndoors.abc  # Ceiling
    abw = PropOpticalIndoors.abw  # Internal wall
    abg = PropOpticalIndoors.abg  # Internal ground
    abm = PropOpticalIndoors.abm  # Internal mass (internal walls)

    # Emissivity [-]
    ec = PropOpticalIndoors.ec  # Ceiling
    eg = PropOpticalIndoors.eg  # Internal wall
    ew = PropOpticalIndoors.ew  # Internal ground
    em = PropOpticalIndoors.em  # Internal mass (internal walls)

    # Heating and cooling parameters
    AC_on = ParHVAC.ACon                         # AC switch (1=on, 0=off)
    Heat_on = ParHVAC.Heatingon                  # Heating switch (1=on, 0=off)
    TspCooling = ParHVAC.TsetpointCooling        # Cooling set-point temperature [K]
    TspHeating = ParHVAC.TsetpointHeating        # Heating set-point temperature [K]
    RHsptCooling = ParHVAC.RHsetpointCooling/100 # Cooling set-point humidity [-]
    RHsptHeating = FT(NaN)                       # Heating set-point humidity (unused)
    COP_AC = ParHVAC.COPAC                       # AC coefficient of performance [-]
    COP_Heat = ParHVAC.COPHeat                   # Heating coefficient of performance [-]
    ACH = ParHVAC.ACH                            # Air changes per hour [1/h]
    f_AC_LEToQ = ParHVAC.f_ACLatentToQ          # Fraction of latent heat to wastewater [0-1]

    AC_onCool = ParHVAC.AC_onCool                # AC cooling switch (1=on, 0=off)
    AC_onDehum = ParHVAC.AC_onDehum              # AC dehumidification switch (1=on, 0=off)

    # Temporal resolution
    dth = ParCalculation.dth  # Time step [h]
    dts = ParCalculation.dts  # Time step [s]

    # Calculate humidity metrics for indoor regulation
    esat_TspCooling =
        FT(611) *
        exp(FT(17.27) * (TspCooling - FT(273.16)) / (FT(237.3) + (TspCooling - FT(273.16))))  # Saturation vapor pressure at set-point [Pa]
    ea_RHspCooling = RHsptCooling * esat_TspCooling  # Vapor pressure at set-point [Pa]
    q_RHspCooling = FT(0.622) * ea_RHspCooling / (Pre - FT(0.378) * ea_RHspCooling)  # Specific humidity at set-point [kg/kg]

    # Calculate constants based on air temperature, pressure, humidity
    L_heat = FT(1000) * (FT(2501.3) - FT(2.361) * (Tatm - FT(273.15)))  # Latent heat [J/kg]
    Cpa = FT(1005) + ((Tatm - FT(273.15) + FT(23.15))^2) / FT(3364)     # Specific heat of air [J/kg K]
    Cpwv = FT(1000) * FT(1.84)                                           # Specific heat of water vapor [J/kg K]
    rho_atm = (Pre / (FT(287.04) * Tatm)) * (1 - (ea/Pre) * (1 - FT(0.622)))  # Dry air density [kg/m³]

    # Building interior energy fluxes
    # Shortwave radiation [W/m²] per surface area
    if ParThermalBuildingInt.IntMassOn
        SWRinB, SWRoutB, SWRabsB = swr_abs_indoors(
            SWRinWsun, SWRinWshd, Hwall, Wroof, abc, abw, abg, abm
        )
    else
        SWRinB, SWRoutB, SWRabsB, SWREBB = swr_abs_indoors_no_int_mass(
            SWRinWsun, SWRinWshd, Hwall, Wroof, abc, abw, abg
        )
        SWRinB = merge(SWRinB, (SWRinInternalMass=zero(FT),))
        SWRoutB = merge(SWRoutB, (SWRoutInternalMass=zero(FT),))
        SWRabsB = merge(SWRabsB, (SWRabsInternalMass=zero(FT),))
    end

    # Longwave radiation [W/m²] per surface area
    # Windows are assumed opaque for longwave radiation
    if ParThermalBuildingInt.IntMassOn
        LWRinB, LWRoutB, LWRabsB = lwr_abs_indoors(
            Tinwallsun,
            Tinwallshd,
            Tceiling,
            Tinground,
            Tintmass,
            Hwall,
            Wroof,
            ec,
            eg,
            em,
            ew,
        )
    else
        LWRinB, LWRoutB, LWRabsB, LWREBB = lwr_abs_indoors_no_int_mass(
            Tinwallsun, Tinwallshd, Tceiling, Tinground, Hwall, Wroof, ec, eg, ew
        )
        LWRinB = merge(LWRinB, (LWRinInternalMass=zero(FT),))
        LWRoutB = merge(LWRoutB, (LWRoutInternalMass=zero(FT),))
        LWRabsB = merge(LWRabsB, (LWRabsInternalMass=zero(FT),))
    end

    # Sensible heat fluxes [W/m²] surface area of each element
    # Positive flux indicates flux from surface to air
    HbinWallSun, HbinWallshd, HBinRoof, HBinGround, HbinIntMass, HbinWindow = sensible_heat_flux_building_interior(
        Tbin, Tinwallsun, Tinwallshd, Tceiling, Tinground, Tintmass, Twindow
    )

    # Calculate total internal sensible heat flux per m² building ground area
    if ParThermalBuildingInt.IntMassOn
        HbuildIn =
            Hwall/Wroof * (HbinWallSun + HbinWallshd) +
            HBinRoof +
            HBinGround +
            2 * Hwall/Wroof * HbinIntMass  # [W/m²] Building ground area
    else
        HbuildIn = Hwall/Wroof * (HbinWallSun + HbinWallshd) + HBinRoof + HBinGround  # [W/m²] Building ground area
        HbinIntMass = zero(FT)
    end

    # Conductive heat flux at building ground/floor [W/m²] ground area
    Gfloor, Tdpfloor = conductive_heat_flux_building_floor(
        Tinground, TingroundDamptm1, Tingroundtm1, ParCalculation, ParThermalBuildingInt
    )

    # Heat storage of internal mass [W/m²] wall height area
    if ParThermalBuildingInt.IntMassOn
        dSinternalMass = heat_storage_change_internal_mass(
            Tintmass, Tintmasstm1, ParThermalBuildingInt, Geometry_m, ParCalculation
        )
    else
        dSinternalMass = zero(FT)
    end

    # Internal sensible and latent heat sources [W/m²] per ground area
    # Positive flux indicates added to indoor air
    Hequip = HVACSchedule.Hequip
    Hpeople = HVACSchedule.Hpeople
    LEequip = HVACSchedule.LEequip
    LEpeople = HVACSchedule.LEpeople

    # Sensible and latent heat load due to ventilation [W/m²] per ground area
    # Positive flux: outdoor air warmer than indoor, heat added from outdoors
    # Negative flux: outdoor air colder than indoor, heat removed from indoors
    Vbuild = Wroof * Hwall
    Hvent = (ACH * dth * Vbuild) / 3600 * Cpa * rho_atm * (Tairout - Tbin) / Wroof
    LEvent = (ACH * dth * Vbuild) / 3600 * rho_atm * L_heat * (qairout - qbin) / Wroof

    # Change in heat and humidity storage in indoor air [W/m²] per ground area
    dSH_air = Vbuild * Cpa * rho_atm * (Tbin - Tbin_tm1) / dts / Wroof
    dSLE_air = Vbuild * rho_atm * L_heat * (qbin - qbin_tm1) / dts / Wroof

    # AC and heating module
    AC_on, AC_onCool, AC_onDehum, Heat_on, H_AC_Heat, LE_AC_Heat = ac_heating_module(
        AC_on,
        Heat_on,
        AC_onCool,
        AC_onDehum,
        ParHVAC,
        HbuildIn,
        Hvent,
        Hequip,
        Hpeople,
        dSH_air,
        LEvent,
        LEequip,
        LEpeople,
        dSLE_air,
    )

    YBuildInt = zeros(FT, 8)  # Initialize energy balance residuals

    # Energy balance for individual surfaces and air volume
    if BEM_on
        # Ceiling energy balance
        YBuildInt[1] = SWRabsB.SWRabsCeiling + LWRabsB.LWRabsCeiling + G2Roof - HBinRoof

        # Sunlit wall energy balance
        YBuildInt[2] =
            SWRabsB.SWRabsWallsun + LWRabsB.LWRabsWallsun + G2WallSun - HbinWallSun

        # Shaded wall energy balance
        YBuildInt[3] =
            SWRabsB.SWRabsWallshd + LWRabsB.LWRabsWallshd + G2WallShade - HbinWallshd

        # Window energy balance (currently not used)
        YBuildInt[4] = TemperatureB[4] - FT(273.15)

        # Ground/floor energy balance
        YBuildInt[5] = SWRabsB.SWRabsGround + LWRabsB.LWRabsGround - Gfloor - HBinGround

        # Internal mass energy balance
        if ParThermalBuildingInt.IntMassOn
            YBuildInt[6] =
                SWRabsB.SWRabsInternalMass + LWRabsB.LWRabsInternalMass - 2*HbinIntMass -
                dSinternalMass
        else
            YBuildInt[6] = Tintmass - FT(273.15)
        end

        # Indoor air temperature and humidity energy balance
        YBuildInt[7] =
            HbuildIn + Hvent + Hequip + Hpeople - dSH_air - H_AC_Heat +
            AC_onCool*FT(1000)*(Tbin - TspCooling) +
            Heat_on*FT(1000)*(Tbin - TspHeating)

        YBuildInt[8] =
            LEvent + LEequip + LEpeople - dSLE_air - LE_AC_Heat +
            AC_onDehum*FT(1e6)*(qbin - q_RHspCooling)
    else
        # If BEM is off, use default values
        YBuildInt[1] = TemperatureB[1] - FT(273.15)
        YBuildInt[2] = TemperatureB[2] - FT(273.15)
        YBuildInt[3] = TemperatureB[3] - FT(273.15)
        YBuildInt[4] = TemperatureB[4] - FT(273.15)
        YBuildInt[5] = TemperatureB[5] - FT(273.15)
        YBuildInt[6] = TemperatureB[6] - FT(273.15)
        YBuildInt[7] = TemperatureB[7] - FT(273.15)
        YBuildInt[8] = FT(1000) * TemperatureB[8] - zero(FT)
    end

    # Waste heat into the canyon, total anthropogenic heat input due to cooling
    # and heating, and building energy demand
    if AC_on
        # Positive if AC is on and added to the canyon as sensible and latent heat source
        # If the moisture removed from the air is fully condensed, we need to
        # account for this additional removal of latent heat in the sensible heat flux
        # The sensible and latent heat flux from the AC accounts for both
        # ventilation heat and internally released heat from walls
        WasteHeat = (;
            SensibleFromAC_Can=(H_AC_Heat + LE_AC_Heat + (H_AC_Heat + LE_AC_Heat)/COP_AC) *
                               Wroof/Wcan,  # [W/m²] canyon ground
            LatentFromAC_Can=zero(FT),  # [W/m²] canyon ground
            WaterFromAC_Can=LE_AC_Heat * Wroof/Wcan,  # [W/m²] canyon ground, water that is condensed and removed as runoff in sewer system
            SensibleFromHeat_Can=zero(FT),  # [W/m²] canyon ground
            LatentFromHeat_Can=zero(FT),  # [W/m²] canyon ground
            SensibleFromVent_Can=-Hvent * Wroof/Wcan,  # [W/m²] canyon ground: negative as hot air is leaving for cooler indoor air
            LatentFromVent_Can=-LEvent * Wroof/Wcan,  # [W/m²] canyon ground
            # Anthropogenic heat added to urban area accounts for internal sources (equipment and people)
            # and additional heat from AC COP not being infinite
            TotAnthInput_URB=(
                (H_AC_Heat + LE_AC_Heat)/COP_AC + Hequip + Hpeople + LEequip + LEpeople
            ) * Wroof/(Wcan + Wroof),  # [W/m²] urban area
        )

        EnergyUse = (;
            EnergyForAC=dth * (H_AC_Heat + LE_AC_Heat) * Wroof/COP_AC,  # [W*h]
            EnergyForAC_H=dth * H_AC_Heat * Wroof/COP_AC,  # [W*h]
            EnergyForAC_LE=dth * LE_AC_Heat * Wroof/COP_AC,  # [W*h]
            EnergyForHeating=zero(FT),  # [W*h]
        )
    elseif Heat_on
        WasteHeat = (
            SensibleFromAC_Can=zero(FT),  # [W/m²] canyon ground
            LatentFromAC_Can=zero(FT),  # [W/m²] canyon ground
            WaterFromAC_Can=zero(FT),  # [W/m²] canyon ground
            SensibleFromHeat_Can=zero(FT),  # [W/m²] canyon ground, all waste heat from heating is added to the canyon through ventilation and conductive heat fluxes through walls
            LatentFromHeat_Can=zero(FT),  # [W/m²] canyon ground
            SensibleFromVent_Can=-Hvent * Wroof/Wcan,  # [W/m²] canyon ground, positive -> heat added to the canyon
            LatentFromVent_Can=-LEvent * Wroof/Wcan,  # [W/m²] canyon ground, positive -> heat added to the canyon
            # Anthropogenic heat added to urban area accounts for internal sources and heating system
            # COP for heating is not considered in waste heat (adding heat is the intended purpose)
            TotAnthInput_URB=(
                -H_AC_Heat - LE_AC_Heat + Hequip + Hpeople + LEequip + LEpeople
            ) * Wroof/(Wcan + Wroof),  # [W/m²] urban area
        )

        EnergyUse = (
            EnergyForAC=zero(FT),  # [W*h]
            EnergyForAC_H=zero(FT),  # [W*h]
            EnergyForAC_LE=zero(FT),  # [W*h]
            EnergyForHeating=dth * (-H_AC_Heat - LE_AC_Heat) * Wroof/COP_Heat,  # [W*h]
        )
    else
        WasteHeat = (
            SensibleFromAC_Can=zero(FT),
            LatentFromAC_Can=zero(FT),
            WaterFromAC_Can=zero(FT),
            SensibleFromHeat_Can=zero(FT),
            LatentFromHeat_Can=zero(FT),
            SensibleFromVent_Can=-Hvent * Wroof/Wcan,
            LatentFromVent_Can=-LEvent * Wroof/Wcan,
            TotAnthInput_URB=(Hequip + Hpeople + LEequip + LEpeople) * Wroof/(Wcan + Wroof),  # [W/m²] urban area
        )

        EnergyUse = (
            EnergyForAC=zero(FT),
            EnergyForAC_H=zero(FT),
            EnergyForAC_LE=zero(FT),
            EnergyForHeating=zero(FT),
        )
    end

    # Apply fraction of Air conditioned rooms for calculation of energy use and waste heat
    WasteHeat = (;
        SensibleFromAC_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.SensibleFromAC_Can),
        LatentFromAC_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.LatentFromAC_Can),
        WaterFromAC_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.WaterFromAC_Can),
        SensibleFromHeat_Can=(
            HVACSchedule.AirConRoomFraction * WasteHeat.SensibleFromHeat_Can
        ),
        LatentFromHeat_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.LatentFromHeat_Can),
        SensibleFromVent_Can=(
            HVACSchedule.AirConRoomFraction * WasteHeat.SensibleFromVent_Can
        ),
        LatentFromVent_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.LatentFromVent_Can),
        TotAnthInput_URB=(HVACSchedule.AirConRoomFraction * WasteHeat.TotAnthInput_URB),
    )

    EnergyUse = (;
        EnergyForAC=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForAC),
        EnergyForAC_H=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForAC_H),
        EnergyForAC_LE=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForAC_LE),
        EnergyForHeating=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForHeating),
    )

    # In case of no BEM

    if BEM_on
        # Prepare energy flux outputs
        # Sensible heat fluxes [W/m²] respective surface area
        HbuildInt = (;
            HBinRoof=HBinRoof,          # interior roof area
            HbinWallSun=HbinWallSun,       # interior sunlit wall area
            HbinWallshd=HbinWallshd,       # interior shaded wall area
            HBinGround=HBinGround,        # interior ground area
            HbinIntMass=HbinIntMass,       # internal mass area (=Wall area)
            HbuildInSurf=HbuildIn,      # building ground area
            Hvent=Hvent,             # building ground area
            Hequip=Hequip,            # building ground area
            Hpeople=Hpeople,           # building ground area
            H_AC_Heat=H_AC_Heat,         # building ground area
            dSH_air=dSH_air,           # building ground area
        )

        # Latent heat fluxes [W/m²] respective surface area
        LEbuildInt = (;
            LEvent=LEvent,            # building ground area
            LEequip=LEequip,           # building ground area
            LEpeople=LEpeople,          # building ground area
            LE_AC_Heat=LE_AC_Heat,        # building ground area
            dSLE_air=dSLE_air,          # building ground area
        )

        # Conductive heat fluxes [W/m²] respective surface area
        GbuildInt = (;
            G2Roof=G2Roof,            # interior roof area
            G2WallSun=G2WallSun,         # interior sunlit wall area
            G2WallShade=G2WallShade,       # interior shaded wall area
            Gfloor=Gfloor,            # interior ground area
            dSinternalMass=dSinternalMass,    # interior wall area
        )

        # Humidity in building interior
        esat_Tbin = 611*exp(17.27*(Tbin-273.15)/(237.3+(Tbin-273.15)))
        e_bin = qbin*Pre/(0.622+0.378*qbin)
        RH_bin = e_bin/esat_Tbin
        HumidityBuilding = (; qbin=qbin, esatbin=esat_Tbin, ebin=e_bin, RHbin=RH_bin)

        # Time varying AC parameters
        ParACHeat = (;
            AC_on=AC_on, AC_onCool=AC_onCool, AC_onDehum=AC_onDehum, Heat_on=Heat_on
        )
    else
        # Prepare energy flux outputs
        # Sensible heat fluxes [W/m²] respective surface area
        HbuildInt = (;
            HBinRoof=zero(FT),          # interior roof area
            HbinWallSun=zero(FT),       # interior sunlit wall area
            HbinWallshd=zero(FT),       # interior shaded wall area
            HBinGround=zero(FT),        # interior ground area
            HbinIntMass=zero(FT),       # internal mass area (=Wall area)
            HbuildInSurf=zero(FT),      # building ground area
            Hvent=zero(FT),             # building ground area
            Hequip=zero(FT),            # building ground area
            Hpeople=zero(FT),           # building ground area
            H_AC_Heat=zero(FT),         # building ground area
            dSH_air=zero(FT),           # building ground area
        )

        # Latent heat fluxes [W/m²] respective surface area
        LEbuildInt = (;
            LEvent=zero(FT),            # building ground area
            LEequip=zero(FT),           # building ground area
            LEpeople=zero(FT),          # building ground area
            LE_AC_Heat=zero(FT),        # building ground area
            dSLE_air=zero(FT),          # building ground area
        )

        # Conductive heat fluxes [W/m²] respective surface area
        GbuildInt = (;
            G2Roof=zero(FT),            # interior roof area
            G2WallSun=zero(FT),         # interior sunlit wall area
            G2WallShade=zero(FT),       # interior shaded wall area
            Gfloor=zero(FT),            # interior ground area
            dSinternalMass=zero(FT),    # interior wall area
        )

        # Humidity in building interior
        esat_Tbin = zero(FT)
        e_bin = zero(FT)
        RH_bin = zero(FT)
        HumidityBuilding = (; qbin=qbin, esatbin=esat_Tbin, ebin=e_bin, RHbin=RH_bin)

        # Time varying AC parameters
        ParACHeat = (; AC_on=false, AC_onCool=false, AC_onDehum=false, Heat_on=false)
    end

    return HbuildInt,
    LEbuildInt,
    GbuildInt,
    SWRabsB,
    LWRabsB,
    Tdpfloor,
    WasteHeat,
    EnergyUse,
    HumidityBuilding,
    ParACHeat,
    YBuildInt
end
function eb_solver_building_output(
    TemperatureC::Vector{FT},
    TemperatureB::Vector{FT},
    TempVecB_ittm::NamedTuple,
    TempVec_ittm::NamedTuple,
    Humidity_ittm::NamedTuple,
    MeteoData::NamedTuple,
    SWRinWsun::FT,
    SWRinWshd::FT,
    G2Roof::FT,
    G2WallSun::FT,
    G2WallShade::FT,
    TempDamp_ittm::NamedTuple,
    SWRabs_t::NamedTuple,
    Geometry_m::ModelComponents.Parameters.UrbanGeometryParameters{FT},
    PropOpticalIndoors::ModelComponents.Parameters.IndoorOpticalProperties{FT},
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    ParCalculation::NamedTuple,
    ParThermalBuildingInt::ModelComponents.Parameters.ThermalBuilding{FT},
    ParWindows::ModelComponents.Parameters.WindowParameters{FT},
    BEM_on::Bool,
    HVACSchedule::NamedTuple,
) where {FT<:AbstractFloat}
    # Initialize building interior temperatures and humidity
    Tceiling = TemperatureB[1]      # Temperature ceiling [K]
    Tinwallsun = TemperatureB[2]    # Temperature sunlit wall [K]
    Tinwallshd = TemperatureB[3]    # Temperature shaded wall [K]
    Twindow = TemperatureB[4]       # Temperature windows [K]
    Tinground = TemperatureB[5]     # Temperature ground [K]
    Tintmass = TemperatureB[6]      # Temperature internal mass [K]
    Tbin = TemperatureB[7]          # Air temperature building interior [K]
    qbin = TemperatureB[8]          # Specific humidity air building interior [kg/kg]

    # External canyon air temperature and specific humidity
    Tairout = TemperatureC[9]       # Air temperature at outdoor calculation height [K]
    qairout = TemperatureC[10]      # Specific humidity at outdoor calculation height [kg/kg]

    # Urban geometry [m]
    Hwall = Geometry_m.Height_canyon  # Building height [m]
    Wroof = Geometry_m.Width_roof     # Roof width [m]
    Wcan = Geometry_m.Width_canyon    # Street/canyon width [m]

    # Temperature [K] and specific humidity [kg/kg] conditions from previous time step
    Tbin_tm1 = TempVecB_ittm.Tbin             # Building internal temperature
    qbin_tm1 = TempVecB_ittm.qbin             # Building internal humidity
    Tingroundtm1 = TempVecB_ittm.Tinground    # Ground temperature
    TingroundDamptm1 = TempDamp_ittm.TDampGroundBuild  # Within ground dampening temperature
    Tintmasstm1 = TempVecB_ittm.Tintmass      # Ground temperature

    Tsurftm1 =
        (
            Wroof * (TempVecB_ittm.Tinground + TempVecB_ittm.Tceiling) +
            Hwall * (TempVecB_ittm.Tinwallshd + TempVecB_ittm.Tinwallsun)
        ) / (2 * Wroof + 2 * Hwall)

    # Atmospheric forcing conditions
    Tatm = MeteoData.Tatm  # Air temperature at atmospheric forcing height [K]
    Pre = MeteoData.Pre    # Air pressure at atmospheric forcing height [Pa]
    ea = MeteoData.ea      # Vapor pressure at atmospheric forcing height [Pa]

    # Albedo [-]
    abc = PropOpticalIndoors.abc  # Ceiling
    abw = PropOpticalIndoors.abw  # Internal wall
    abg = PropOpticalIndoors.abg  # Internal ground
    abm = PropOpticalIndoors.abm  # Internal mass (internal walls)

    # Emissivity [-]
    ec = PropOpticalIndoors.ec  # Ceiling
    eg = PropOpticalIndoors.eg  # Internal wall
    ew = PropOpticalIndoors.ew  # Internal ground
    em = PropOpticalIndoors.em  # Internal mass (internal walls)

    # Heating and cooling parameters
    AC_on = ParHVAC.ACon                         # AC switch (1=on, 0=off)
    Heat_on = ParHVAC.Heatingon                  # Heating switch (1=on, 0=off)
    TspCooling = ParHVAC.TsetpointCooling        # Cooling set-point temperature [K]
    TspHeating = ParHVAC.TsetpointHeating        # Heating set-point temperature [K]
    RHsptCooling = ParHVAC.RHsetpointCooling/100 # Cooling set-point humidity [-]
    RHsptHeating = FT(NaN)                       # Heating set-point humidity (unused)
    COP_AC = ParHVAC.COPAC                       # AC coefficient of performance [-]
    COP_Heat = ParHVAC.COPHeat                   # Heating coefficient of performance [-]
    ACH = ParHVAC.ACH                            # Air changes per hour [1/h]
    f_AC_LEToQ = ParHVAC.f_ACLatentToQ          # Fraction of latent heat to wastewater [0-1]

    AC_onCool = ParHVAC.AC_onCool                # AC cooling switch (1=on, 0=off)
    AC_onDehum = ParHVAC.AC_onDehum              # AC dehumidification switch (1=on, 0=off)

    # Temporal resolution
    dth = ParCalculation.dth  # Time step [h]
    dts = ParCalculation.dts  # Time step [s]

    # Calculate humidity metrics for indoor regulation
    esat_TspCooling =
        FT(611) *
        exp(FT(17.27) * (TspCooling - FT(273.16)) / (FT(237.3) + (TspCooling - FT(273.16))))  # Saturation vapor pressure at set-point [Pa]
    ea_RHspCooling = RHsptCooling * esat_TspCooling  # Vapor pressure at set-point [Pa]
    q_RHspCooling = FT(0.622) * ea_RHspCooling / (Pre - FT(0.378) * ea_RHspCooling)  # Specific humidity at set-point [kg/kg]

    # Calculate constants based on air temperature, pressure, humidity
    L_heat = FT(1000) * (FT(2501.3) - FT(2.361) * (Tatm - FT(273.15)))  # Latent heat [J/kg]
    Cpa = FT(1005) + ((Tatm - FT(273.15) + FT(23.15))^2) / FT(3364)     # Specific heat of air [J/kg K]
    Cpwv = FT(1000) * FT(1.84)                                           # Specific heat of water vapor [J/kg K]
    rho_atm = (Pre / (FT(287.04) * Tatm)) * (1 - (ea/Pre) * (1 - FT(0.622)))  # Dry air density [kg/m³]

    # Building interior energy fluxes
    # Shortwave radiation [W/m²] per surface area
    if ParThermalBuildingInt.IntMassOn
        SWRinB, SWRoutB, SWRabsB = swr_abs_indoors(
            SWRinWsun, SWRinWshd, Hwall, Wroof, abc, abw, abg, abm
        )
    else
        SWRinB, SWRoutB, SWRabsB, SWREBB = swr_abs_indoors_no_int_mass(
            SWRinWsun, SWRinWshd, Hwall, Wroof, abc, abw, abg
        )
        SWRinB = merge(SWRinB, (SWRinInternalMass=zero(FT),))
        SWRoutB = merge(SWRoutB, (SWRoutInternalMass=zero(FT),))
        SWRabsB = merge(SWRabsB, (SWRabsInternalMass=zero(FT),))
    end

    # Longwave radiation [W/m²] per surface area
    # Windows are assumed opaque for longwave radiation
    if ParThermalBuildingInt.IntMassOn
        LWRinB, LWRoutB, LWRabsB = lwr_abs_indoors(
            Tinwallsun,
            Tinwallshd,
            Tceiling,
            Tinground,
            Tintmass,
            Hwall,
            Wroof,
            ec,
            eg,
            em,
            ew,
        )
    else
        LWRinB, LWRoutB, LWRabsB, LWREBB = lwr_abs_indoors_no_int_mass(
            Tinwallsun, Tinwallshd, Tceiling, Tinground, Hwall, Wroof, ec, eg, ew
        )
        LWRinB = merge(LWRinB, (LWRinInternalMass=zero(FT),))
        LWRoutB = merge(LWRoutB, (LWRoutInternalMass=zero(FT),))
        LWRabsB = merge(LWRabsB, (LWRabsInternalMass=zero(FT),))
    end

    # Sensible heat fluxes [W/m²] surface area of each element
    # Positive flux indicates flux from surface to air
    HbinWallSun, HbinWallshd, HBinRoof, HBinGround, HbinIntMass, HbinWindow = sensible_heat_flux_building_interior(
        Tbin, Tinwallsun, Tinwallshd, Tceiling, Tinground, Tintmass, Twindow
    )

    # Calculate total internal sensible heat flux per m² building ground area
    if ParThermalBuildingInt.IntMassOn
        HbuildIn =
            Hwall/Wroof * (HbinWallSun + HbinWallshd) +
            HBinRoof +
            HBinGround +
            2 * Hwall/Wroof * HbinIntMass  # [W/m²] Building ground area
    else
        HbuildIn = Hwall/Wroof * (HbinWallSun + HbinWallshd) + HBinRoof + HBinGround  # [W/m²] Building ground area
        HbinIntMass = zero(FT)
    end

    # Conductive heat flux at building ground/floor [W/m²] ground area
    Gfloor, Tdpfloor = conductive_heat_flux_building_floor(
        Tinground, TingroundDamptm1, Tingroundtm1, ParCalculation, ParThermalBuildingInt
    )

    # Heat storage of internal mass [W/m²] wall height area
    if ParThermalBuildingInt.IntMassOn
        dSinternalMass = heat_storage_change_internal_mass(
            Tintmass, Tintmasstm1, ParThermalBuildingInt, Geometry_m, ParCalculation
        )
    else
        dSinternalMass = zero(FT)
    end

    # Internal sensible and latent heat sources [W/m²] per ground area
    # Positive flux indicates added to indoor air
    Hequip = HVACSchedule.Hequip
    Hpeople = HVACSchedule.Hpeople
    LEequip = HVACSchedule.LEequip
    LEpeople = HVACSchedule.LEpeople

    # Sensible and latent heat load due to ventilation [W/m²] per ground area
    # Positive flux: outdoor air warmer than indoor, heat added from outdoors
    # Negative flux: outdoor air colder than indoor, heat removed from indoors
    Vbuild = Wroof * Hwall
    Hvent = (ACH * dth * Vbuild) / 3600 * Cpa * rho_atm * (Tairout - Tbin) / Wroof
    LEvent = (ACH * dth * Vbuild) / 3600 * rho_atm * L_heat * (qairout - qbin) / Wroof

    # Change in heat and humidity storage in indoor air [W/m²] per ground area
    dSH_air = Vbuild * Cpa * rho_atm * (Tbin - Tbin_tm1) / dts / Wroof
    dSLE_air = Vbuild * rho_atm * L_heat * (qbin - qbin_tm1) / dts / Wroof

    # AC and heating module
    AC_on, AC_onCool, AC_onDehum, Heat_on, H_AC_Heat, LE_AC_Heat = ac_heating_module(
        AC_on,
        Heat_on,
        AC_onCool,
        AC_onDehum,
        ParHVAC,
        HbuildIn,
        Hvent,
        Hequip,
        Hpeople,
        dSH_air,
        LEvent,
        LEequip,
        LEpeople,
        dSLE_air,
    )

    YBuildInt = zeros(FT, 8)  # Initialize energy balance residuals

    # Energy balance for individual surfaces and air volume
    if BEM_on
        # Ceiling energy balance
        YBuildInt[1] = SWRabsB.SWRabsCeiling + LWRabsB.LWRabsCeiling + G2Roof - HBinRoof

        # Sunlit wall energy balance
        YBuildInt[2] =
            SWRabsB.SWRabsWallsun + LWRabsB.LWRabsWallsun + G2WallSun - HbinWallSun

        # Shaded wall energy balance
        YBuildInt[3] =
            SWRabsB.SWRabsWallshd + LWRabsB.LWRabsWallshd + G2WallShade - HbinWallshd

        # Window energy balance (currently not used)
        YBuildInt[4] = TemperatureB[4] - FT(273.15)

        # Ground/floor energy balance
        YBuildInt[5] = SWRabsB.SWRabsGround + LWRabsB.LWRabsGround - Gfloor - HBinGround

        # Internal mass energy balance
        if ParThermalBuildingInt.IntMassOn
            YBuildInt[6] =
                SWRabsB.SWRabsInternalMass + LWRabsB.LWRabsInternalMass - 2*HbinIntMass -
                dSinternalMass
        else
            YBuildInt[6] = Tintmass - FT(273.15)
        end

        # Indoor air temperature and humidity energy balance
        YBuildInt[7] =
            HbuildIn + Hvent + Hequip + Hpeople - dSH_air - H_AC_Heat +
            AC_onCool*FT(1000)*(Tbin - TspCooling) +
            Heat_on*FT(1000)*(Tbin - TspHeating)

        YBuildInt[8] =
            LEvent + LEequip + LEpeople - dSLE_air - LE_AC_Heat +
            AC_onDehum*FT(1e6)*(qbin - q_RHspCooling)
    else
        # If BEM is off, use default values
        YBuildInt[1] = TemperatureB[1] - FT(273.15)
        YBuildInt[2] = TemperatureB[2] - FT(273.15)
        YBuildInt[3] = TemperatureB[3] - FT(273.15)
        YBuildInt[4] = TemperatureB[4] - FT(273.15)
        YBuildInt[5] = TemperatureB[5] - FT(273.15)
        YBuildInt[6] = TemperatureB[6] - FT(273.15)
        YBuildInt[7] = TemperatureB[7] - FT(273.15)
        YBuildInt[8] = FT(1000) * TemperatureB[8] - zero(FT)
    end

    # Waste heat into the canyon, total anthropogenic heat input due to cooling
    # and heating, and building energy demand
    if AC_on
        # Positive if AC is on and added to the canyon as sensible and latent heat source
        # If the moisture removed from the air is fully condensed, we need to
        # account for this additional removal of latent heat in the sensible heat flux
        # The sensible and latent heat flux from the AC accounts for both
        # ventilation heat and internally released heat from walls
        WasteHeat = (;
            SensibleFromAC_Can=(H_AC_Heat + LE_AC_Heat + (H_AC_Heat + LE_AC_Heat)/COP_AC) *
                               Wroof/Wcan,  # [W/m²] canyon ground
            LatentFromAC_Can=zero(FT),  # [W/m²] canyon ground
            WaterFromAC_Can=LE_AC_Heat * Wroof/Wcan,  # [W/m²] canyon ground, water that is condensed and removed as runoff in sewer system
            SensibleFromHeat_Can=zero(FT),  # [W/m²] canyon ground
            LatentFromHeat_Can=zero(FT),  # [W/m²] canyon ground
            SensibleFromVent_Can=-Hvent * Wroof/Wcan,  # [W/m²] canyon ground: negative as hot air is leaving for cooler indoor air
            LatentFromVent_Can=-LEvent * Wroof/Wcan,  # [W/m²] canyon ground
            # Anthropogenic heat added to urban area accounts for internal sources (equipment and people)
            # and additional heat from AC COP not being infinite
            TotAnthInput_URB=(
                (H_AC_Heat + LE_AC_Heat)/COP_AC + Hequip + Hpeople + LEequip + LEpeople
            ) * Wroof/(Wcan + Wroof),  # [W/m²] urban area
        )

        EnergyUse = (;
            EnergyForAC=dth * (H_AC_Heat + LE_AC_Heat) * Wroof/COP_AC,  # [W*h]
            EnergyForAC_H=dth * H_AC_Heat * Wroof/COP_AC,  # [W*h]
            EnergyForAC_LE=dth * LE_AC_Heat * Wroof/COP_AC,  # [W*h]
            EnergyForHeating=zero(FT),  # [W*h]
        )
    elseif Heat_on
        WasteHeat = (
            SensibleFromAC_Can=zero(FT),  # [W/m²] canyon ground
            LatentFromAC_Can=zero(FT),  # [W/m²] canyon ground
            WaterFromAC_Can=zero(FT),  # [W/m²] canyon ground
            SensibleFromHeat_Can=zero(FT),  # [W/m²] canyon ground, all waste heat from heating is added to the canyon through ventilation and conductive heat fluxes through walls
            LatentFromHeat_Can=zero(FT),  # [W/m²] canyon ground
            SensibleFromVent_Can=-Hvent * Wroof/Wcan,  # [W/m²] canyon ground, positive -> heat added to the canyon
            LatentFromVent_Can=-LEvent * Wroof/Wcan,  # [W/m²] canyon ground, positive -> heat added to the canyon
            # Anthropogenic heat added to urban area accounts for internal sources and heating system
            # COP for heating is not considered in waste heat (adding heat is the intended purpose)
            TotAnthInput_URB=(
                -H_AC_Heat - LE_AC_Heat + Hequip + Hpeople + LEequip + LEpeople
            ) * Wroof/(Wcan + Wroof),  # [W/m²] urban area
        )

        EnergyUse = (
            EnergyForAC=zero(FT),  # [W*h]
            EnergyForAC_H=zero(FT),  # [W*h]
            EnergyForAC_LE=zero(FT),  # [W*h]
            EnergyForHeating=dth * (-H_AC_Heat - LE_AC_Heat) * Wroof/COP_Heat,  # [W*h]
        )
    else
        WasteHeat = (
            SensibleFromAC_Can=zero(FT),
            LatentFromAC_Can=zero(FT),
            WaterFromAC_Can=zero(FT),
            SensibleFromHeat_Can=zero(FT),
            LatentFromHeat_Can=zero(FT),
            SensibleFromVent_Can=-Hvent * Wroof/Wcan,
            LatentFromVent_Can=-LEvent * Wroof/Wcan,
            TotAnthInput_URB=(Hequip + Hpeople + LEequip + LEpeople) * Wroof/(Wcan + Wroof),  # [W/m²] urban area
        )

        EnergyUse = (
            EnergyForAC=zero(FT),
            EnergyForAC_H=zero(FT),
            EnergyForAC_LE=zero(FT),
            EnergyForHeating=zero(FT),
        )
    end

    # Apply fraction of Air conditioned rooms for calculation of energy use and waste heat
    WasteHeat = (;
        SensibleFromAC_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.SensibleFromAC_Can),
        LatentFromAC_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.LatentFromAC_Can),
        WaterFromAC_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.WaterFromAC_Can),
        SensibleFromHeat_Can=(
            HVACSchedule.AirConRoomFraction * WasteHeat.SensibleFromHeat_Can
        ),
        LatentFromHeat_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.LatentFromHeat_Can),
        SensibleFromVent_Can=(
            HVACSchedule.AirConRoomFraction * WasteHeat.SensibleFromVent_Can
        ),
        LatentFromVent_Can=(HVACSchedule.AirConRoomFraction * WasteHeat.LatentFromVent_Can),
        TotAnthInput_URB=(HVACSchedule.AirConRoomFraction * WasteHeat.TotAnthInput_URB),
    )

    EnergyUse = (;
        EnergyForAC=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForAC),
        EnergyForAC_H=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForAC_H),
        EnergyForAC_LE=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForAC_LE),
        EnergyForHeating=(HVACSchedule.AirConRoomFraction * EnergyUse.EnergyForHeating),
    )

    # In case of no BEM

    if BEM_on
        # Prepare energy flux outputs
        # Sensible heat fluxes [W/m²] respective surface area
        HbuildInt = (;
            HBinRoof=HBinRoof,          # interior roof area
            HbinWallSun=HbinWallSun,       # interior sunlit wall area
            HbinWallshd=HbinWallshd,       # interior shaded wall area
            HBinGround=HBinGround,        # interior ground area
            HbinIntMass=HbinIntMass,       # internal mass area (=Wall area)
            HbuildInSurf=HbuildIn,      # building ground area
            Hvent=Hvent,             # building ground area
            Hequip=Hequip,            # building ground area
            Hpeople=Hpeople,           # building ground area
            H_AC_Heat=H_AC_Heat,         # building ground area
            dSH_air=dSH_air,           # building ground area
        )

        # Latent heat fluxes [W/m²] respective surface area
        LEbuildInt = (;
            LEvent=LEvent,            # building ground area
            LEequip=LEequip,           # building ground area
            LEpeople=LEpeople,          # building ground area
            LE_AC_Heat=LE_AC_Heat,        # building ground area
            dSLE_air=dSLE_air,          # building ground area
        )

        # Conductive heat fluxes [W/m²] respective surface area
        GbuildInt = (;
            G2Roof=G2Roof,            # interior roof area
            G2WallSun=G2WallSun,         # interior sunlit wall area
            G2WallShade=G2WallShade,       # interior shaded wall area
            Gfloor=Gfloor,            # interior ground area
            dSinternalMass=dSinternalMass,    # interior wall area
        )

        # Humidity in building interior
        esat_Tbin = 611*exp(17.27*(Tbin-273.15)/(237.3+(Tbin-273.15)))
        e_bin = qbin*Pre/(0.622+0.378*qbin)
        RH_bin = e_bin/esat_Tbin
        HumidityBuilding = (; qbin=qbin, esatbin=esat_Tbin, ebin=e_bin, RHbin=RH_bin)

        # Time varying AC parameters
        ParACHeat = (;
            AC_on=AC_on, AC_onCool=AC_onCool, AC_onDehum=AC_onDehum, Heat_on=Heat_on
        )
    else
        # Prepare energy flux outputs
        # Sensible heat fluxes [W/m²] respective surface area
        HbuildInt = (;
            HBinRoof=zero(FT),          # interior roof area
            HbinWallSun=zero(FT),       # interior sunlit wall area
            HbinWallshd=zero(FT),       # interior shaded wall area
            HBinGround=zero(FT),        # interior ground area
            HbinIntMass=zero(FT),       # internal mass area (=Wall area)
            HbuildInSurf=zero(FT),      # building ground area
            Hvent=zero(FT),             # building ground area
            Hequip=zero(FT),            # building ground area
            Hpeople=zero(FT),           # building ground area
            H_AC_Heat=zero(FT),         # building ground area
            dSH_air=zero(FT),           # building ground area
        )

        # Latent heat fluxes [W/m²] respective surface area
        LEbuildInt = (;
            LEvent=zero(FT),            # building ground area
            LEequip=zero(FT),           # building ground area
            LEpeople=zero(FT),          # building ground area
            LE_AC_Heat=zero(FT),        # building ground area
            dSLE_air=zero(FT),          # building ground area
        )

        # Conductive heat fluxes [W/m²] respective surface area
        GbuildInt = (;
            G2Roof=zero(FT),            # interior roof area
            G2WallSun=zero(FT),         # interior sunlit wall area
            G2WallShade=zero(FT),       # interior shaded wall area
            Gfloor=zero(FT),            # interior ground area
            dSinternalMass=zero(FT),    # interior wall area
        )

        # Humidity in building interior
        esat_Tbin = zero(FT)
        e_bin = zero(FT)
        RH_bin = zero(FT)
        HumidityBuilding = (; qbin=qbin, esatbin=esat_Tbin, ebin=e_bin, RHbin=RH_bin)

        # Time varying AC parameters
        ParACHeat = (; AC_on=false, AC_onCool=false, AC_onDehum=false, Heat_on=false)
    end

    return HbuildInt,
    LEbuildInt,
    GbuildInt,
    SWRabsB,
    LWRabsB,
    Tdpfloor,
    WasteHeat,
    EnergyUse,
    HumidityBuilding,
    ParACHeat,
    YBuildInt
end
