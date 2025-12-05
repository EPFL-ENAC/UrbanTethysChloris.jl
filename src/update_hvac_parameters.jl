function update_hvac_parameters(
    ParHVACorig::ModelComponents.Parameters.HVACParameters{FT},
    ParHVAC::ModelComponents.Parameters.HVACParameters{FT},
    EnergyUse::NamedTuple,
    Tbin::FT,
    qbin::FT,
) where {FT<:AbstractFloat}
    ACon = ParHVAC.ACon
    AC_onCool = ParHVAC.AC_onCool
    AC_onDehum = ParHVAC.AC_onDehum
    Heatingon = ParHVAC.Heatingon
    MasterOn = ParHVAC.MasterOn

    Tbin_r = round(Tbin; digits=4)
    qbin_r = round(qbin; digits=8)

    if ParHVACorig.ACon && Tbin_r>(ParHVAC.TsetpointCooling+0.01) ||
        qbin_r>(ParHVAC.q_RHspCooling+10^-6)
        ACon = 1;
        AC_onCool = 1;
        AC_onDehum = 1;
        Heatingon = 0;
        MasterOn = 1;
        if qbin_r<(ParHVAC.q_RHspCooling+10^-6) && round(EnergyUse.EnergyForAC_LE, 1)==0
            AC_onDehum = 0;
        elseif Tbin_r<(ParHVAC.TsetpointCooling+0.01) &&
            round(EnergyUse.EnergyForAC_H, 1)==0
            AC_onCool = 0;
        end
    elseif ParHVACorig.ACon &&
        EnergyUse.EnergyForAC_H>0 &&
        qbin_r>(ParHVAC.q_RHspCooling+10^-6)
        ACon = 1;
        AC_onCool = 1;
        AC_onDehum = 1;
        MasterOn = 1;
    elseif ParHVACorig.Heatingon &&
        Tbin_r<(ParHVAC.TsetpointHeating-0.01) &&
        round(EnergyUse.EnergyForHeating)==0
        ACon = 0;
        AC_onCool = 0;
        AC_onDehum = 0;
        Heatingon = 1;
        MasterOn = 1;
    end

    if ParHVACorig.ACon==1 && EnergyUse.EnergyForAC_H<-10^-6 ||
        EnergyUse.EnergyForAC_LE<-10^-6
        if EnergyUse.EnergyForAC_H<-10^-6 && EnergyUse.EnergyForAC_LE<-10^-6
            ACon = 0;
            AC_onCool = 0;
            AC_onDehum = 0;
            MasterOn = 1;
        elseif EnergyUse.EnergyForAC_H<-10^-6 && EnergyUse.EnergyForAC_LE>10^-6
            ACon = 1;
            AC_onCool = 0;
            AC_onDehum = 1;
            MasterOn = 1;
        elseif EnergyUse.EnergyForAC_H>10^-6 && EnergyUse.EnergyForAC_LE<-10^-6
            ACon = 1;
            AC_onCool = 1;
            AC_onDehum = 0;
            MasterOn = 1;
        elseif EnergyUse.EnergyForAC_H<-10^-6 && EnergyUse.EnergyForAC_LE==0
            ACon = 0;
            AC_onCool = 0;
            AC_onDehum = 0;
            MasterOn = 1;
        elseif EnergyUse.EnergyForAC_LE<-10^-6 && EnergyUse.EnergyForAC_H==0
            ACon = 0;
            AC_onCool = 0;
            AC_onDehum = 0;
            MasterOn = 1;
        end
    elseif ParHVACorig.Heatingon==1 && EnergyUse.EnergyForHeating<-10^-6
        Heatingon = 0;
        MasterOn = 1;
    end

    if ACon ≠ ParHVAC.ACon
        @info "ACon switched from $(ParHVAC.ACon) to $ACon"
    end

    if AC_onCool ≠ ParHVAC.AC_onCool
        @info "AC_onCool switched from $(ParHVAC.ACon) to $ACon"
    end

    if AC_onDehum ≠ ParHVAC.AC_onDehum
        @info "AC_onDehum switched from $(ParHVAC.AC_onDehum) to $AC_onDehum"
    end

    if Heatingon ≠ ParHVAC.Heatingon
        @info "Heatingon switched from $(ParHVAC.Heatingon) to $Heatingon"
    end

    if MasterOn ≠ ParHVAC.MasterOn
        @info "MasterOn switched from $(ParHVAC.ACon) to $ACon"
    end

    return ModelComponents.Parameters.HVACParameters{FT}(;
        ACon=ACon,
        AC_onCool=AC_onCool,
        AC_onDehum=AC_onDehum,
        MasterOn=MasterOn,
        Heatingon=Heatingon,
        TsetpointCooling=ParHVAC.TsetpointCooling,
        TsetpointHeating=ParHVAC.TsetpointHeating,
        RHsetpointCooling=ParHVAC.RHsetpointCooling,
        RHsetpointHeating=ParHVAC.RHsetpointHeating,
        ACH=ParHVAC.ACH,
        COPAC=ParHVAC.COPAC,
        COPHeat=ParHVAC.COPHeat,
        f_ACLatentToQ=ParHVAC.f_ACLatentToQ,
        q_RHspCooling=ParHVAC.q_RHspCooling,
    )
end
