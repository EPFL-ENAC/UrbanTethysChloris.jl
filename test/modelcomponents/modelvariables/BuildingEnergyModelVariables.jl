using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    TempVecB,
    initialize_tempvecb,
    HumidityBuilding,
    initialize_humidity_building,
    HbuildInt,
    initialize_hbuildint,
    LEbuildInt,
    initialize_lebuildint,
    GbuildInt,
    initialize_gbuildint,
    SWRabsB,
    initialize_swrabsb,
    LWRabsB,
    initialize_lwrabsb,
    BEMWasteHeat,
    initialize_bemwasteheat,
    BEMEnergyUse,
    initialize_bemenergyuse,
    ParACHeat_ts,
    initialize_paracheat_ts,
    BuildingEnergyModelVariables,
    initialize_building_energy_model_variables,
    TimeSlice,
    TimeSeries

FT = Float64

@testset "Subsets" begin
    @testset "TempVecB" begin
        @testset "TimeSlice" begin
            tempvecb = initialize_tempvecb(FT, TimeSlice())

            # Test structure
            @test tempvecb isa TempVecB{FT,0}

            # Test all fields are accessible
            for field in fieldnames(TempVecB)
                @test isa(getproperty(tempvecb, field), FT)
                @test getproperty(tempvecb, field) == 0.0
            end
        end

        @testset "TimeSeries" begin
            Tatm = FT(300.0)  # Initial temperature
            AtmSpecific = FT(0.01)  # Initial specific humidity

            hours = 24
            tempvecb = initialize_tempvecb(FT, TimeSeries(), hours, Tatm, AtmSpecific)

            @test tempvecb isa TempVecB{FT,1}

            # Test initial values
            @test tempvecb.Tceiling[1] == Tatm
            @test tempvecb.Tinwallsun[1] == Tatm
            @test tempvecb.qbin[1] == AtmSpecific

            # Test all fields are accessible and have correct dimensions
            for field in fieldnames(TempVecB)
                @test isa(getproperty(tempvecb, field), Array{FT,1})
                @test size(getproperty(tempvecb, field)) == (hours,)
            end
        end
    end

    @testset "HumidityBuilding" begin
        @testset "TimeSlice" begin
            humidityb = initialize_humidity_building(FT, TimeSlice())

            # Test structure
            @test humidityb isa HumidityBuilding{FT,0}

            # Test field access for scalar case
            @test humidityb.esatbin === 0.0
            @test humidityb.ebin === 0.0
            @test humidityb.RHbin === 0.0

            # Test all fields are accessible
            for field in fieldnames(HumidityBuilding)
                @test isa(getproperty(humidityb, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            humidityb = initialize_humidity_building(FT, TimeSeries(), hours)

            @test humidityb isa HumidityBuilding{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(HumidityBuilding)
                @test isa(getproperty(humidityb, field), Array{FT,1})
                @test size(getproperty(humidityb, field)) == (hours,)
                @test getproperty(humidityb, field)[1] == 0
            end
        end
    end

    @testset "HbuildInt" begin
        @testset "TimeSlice" begin
            hbuildint = initialize_hbuildint(FT, TimeSlice())

            # Test structure
            @test hbuildint isa HbuildInt{FT,0}

            # Test field access for scalar case
            @test hbuildint.HBinRoof === 0.0
            @test hbuildint.HbuildInSurf === 0.0
            @test hbuildint.dSH_air === 0.0

            # Test all fields are accessible
            for field in fieldnames(HbuildInt)
                @test isa(getproperty(hbuildint, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            hbuildint = initialize_hbuildint(FT, TimeSeries(), hours)

            @test hbuildint isa HbuildInt{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(HbuildInt)
                @test isa(getproperty(hbuildint, field), Array{FT,1})
                @test size(getproperty(hbuildint, field)) == (hours,)
                @test getproperty(hbuildint, field)[1] == 0
            end
        end
    end

    @testset "LEbuildInt" begin
        @testset "TimeSlice" begin
            lebuildint = initialize_lebuildint(FT, TimeSlice())

            # Test structure
            @test lebuildint isa LEbuildInt{FT,0}

            # Test field access for scalar case
            @test lebuildint.LEvent === 0.0
            @test lebuildint.LEpeople === 0.0
            @test lebuildint.dSLE_air === 0.0

            # Test all fields are accessible
            for field in fieldnames(LEbuildInt)
                @test isa(getproperty(lebuildint, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            lebuildint = initialize_lebuildint(FT, TimeSeries(), hours)

            @test lebuildint isa LEbuildInt{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(LEbuildInt)
                @test isa(getproperty(lebuildint, field), Array{FT,1})
                @test size(getproperty(lebuildint, field)) == (hours,)
                @test getproperty(lebuildint, field)[1] == 0
            end
        end
    end

    @testset "GbuildInt" begin
        @testset "TimeSlice" begin
            gbuildint = initialize_gbuildint(FT, TimeSlice())

            # Test structure
            @test gbuildint isa GbuildInt{FT,0}

            # Test field access for scalar case
            @test gbuildint.G2Roof === 0.0
            @test gbuildint.G2WallSun === 0.0
            @test gbuildint.dSinternalMass === 0.0

            # Test all fields are accessible
            for field in fieldnames(GbuildInt)
                @test isa(getproperty(gbuildint, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            gbuildint = initialize_gbuildint(FT, TimeSeries(), hours)

            @test gbuildint isa GbuildInt{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(GbuildInt)
                @test isa(getproperty(gbuildint, field), Array{FT,1})
                @test size(getproperty(gbuildint, field)) == (hours,)
                @test getproperty(gbuildint, field)[1] == 0
            end
        end
    end

    @testset "SWRabsB" begin
        @testset "TimeSlice" begin
            swrabsb = initialize_swrabsb(FT, TimeSlice())

            # Test structure
            @test swrabsb isa SWRabsB{FT,0}

            # Test field access for scalar case
            @test swrabsb.SWRabsCeiling === 0.0
            @test swrabsb.SWRabsWallsun === 0.0
            @test swrabsb.SWRabsInternalMass === 0.0

            # Test all fields are accessible
            for field in fieldnames(SWRabsB)
                @test isa(getproperty(swrabsb, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            swrabsb = initialize_swrabsb(FT, TimeSeries(), hours)

            @test swrabsb isa SWRabsB{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(SWRabsB)
                @test isa(getproperty(swrabsb, field), Array{FT,1})
                @test size(getproperty(swrabsb, field)) == (hours,)
                @test getproperty(swrabsb, field)[1] == 0
            end
        end
    end

    @testset "LWRabsB" begin
        @testset "TimeSlice" begin
            lwrabsb = initialize_lwrabsb(FT, TimeSlice())

            # Test structure
            @test lwrabsb isa LWRabsB{FT,0}

            # Test field access for scalar case
            @test lwrabsb.LWRabsCeiling === 0.0
            @test lwrabsb.LWRabsWallsun === 0.0
            @test lwrabsb.LWRabsInternalMass === 0.0

            # Test all fields are accessible
            for field in fieldnames(LWRabsB)
                @test isa(getproperty(lwrabsb, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            lwrabsb = initialize_lwrabsb(FT, TimeSeries(), hours)

            @test lwrabsb isa LWRabsB{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(LWRabsB)
                @test isa(getproperty(lwrabsb, field), Array{FT,1})
                @test size(getproperty(lwrabsb, field)) == (hours,)
                @test getproperty(lwrabsb, field)[1] == 0
            end
        end
    end

    @testset "BEMWasteHeat" begin
        @testset "TimeSlice" begin
            bemwasteheat = initialize_bemwasteheat(FT, TimeSlice())

            # Test structure
            @test bemwasteheat isa BEMWasteHeat{FT,0}

            # Test field access for scalar case
            @test bemwasteheat.SensibleFromAC_Can === 0.0
            @test bemwasteheat.LatentFromHeat_Can === 0.0
            @test bemwasteheat.TotAnthInput_URB === 0.0

            # Test all fields are accessible
            for field in fieldnames(BEMWasteHeat)
                @test isa(getproperty(bemwasteheat, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            bemwasteheat = initialize_bemwasteheat(FT, TimeSeries(), hours)

            @test bemwasteheat isa BEMWasteHeat{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(BEMWasteHeat)
                @test isa(getproperty(bemwasteheat, field), Array{FT,1})
                @test size(getproperty(bemwasteheat, field)) == (hours,)
                @test getproperty(bemwasteheat, field)[1] == 0
            end
        end
    end

    @testset "BEMEnergyUse" begin
        @testset "TimeSlice" begin
            bemenergyuse = initialize_bemenergyuse(FT, TimeSlice())

            # Test structure
            @test bemenergyuse isa BEMEnergyUse{FT,0}

            # Test field access for scalar case
            @test bemenergyuse.EnergyForAC === 0.0
            @test bemenergyuse.EnergyForAC_LE === 0.0
            @test bemenergyuse.EnergyForHeating === 0.0

            # Test all fields are accessible
            for field in fieldnames(BEMEnergyUse)
                @test isa(getproperty(bemenergyuse, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            bemenergyuse = initialize_bemenergyuse(FT, TimeSeries(), hours)

            @test bemenergyuse isa BEMEnergyUse{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(BEMEnergyUse)
                @test isa(getproperty(bemenergyuse, field), Array{FT,1})
                @test size(getproperty(bemenergyuse, field)) == (hours,)
                @test getproperty(bemenergyuse, field)[1] == 0
            end
        end
    end

    @testset "ParACHeat_ts" begin
        @testset "TimeSlice" begin
            paracheat_ts = initialize_paracheat_ts(FT, TimeSlice())

            # Test structure
            @test paracheat_ts isa ParACHeat_ts{FT,0}

            # Test field access for scalar case
            @test paracheat_ts.AC_on === 0.0
            @test paracheat_ts.AC_onCool === 0.0
            @test paracheat_ts.Heat_on === 0.0

            # Test all fields are accessible
            for field in fieldnames(ParACHeat_ts)
                @test isa(getproperty(paracheat_ts, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            paracheat_ts = initialize_paracheat_ts(FT, TimeSeries(), hours)

            @test paracheat_ts isa ParACHeat_ts{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(ParACHeat_ts)
                @test isa(getproperty(paracheat_ts, field), Array{FT,1})
                @test size(getproperty(paracheat_ts, field)) == (hours,)
                @test getproperty(paracheat_ts, field)[1] == 0
            end
        end
    end
end

@testset "BuildingEnergyModelVariables" begin
    @testset "TimeSlice" begin
        bem_vars = initialize_building_energy_model_variables(FT, TimeSlice())

        # Test structure
        @test bem_vars isa BuildingEnergyModelVariables{FT,0}

        # Test fields initialized to zero
        @test isa(bem_vars.TempVecB, TempVecB{FT,0})
        @test isa(bem_vars.HumidityBuilding, HumidityBuilding{FT,0})
        @test isa(bem_vars.HbuildInt, HbuildInt{FT,0})
        @test isa(bem_vars.LEbuildInt, LEbuildInt{FT,0})
        @test isa(bem_vars.GbuildInt, GbuildInt{FT,0})
        @test isa(bem_vars.SWRabsB, SWRabsB{FT,0})
        @test isa(bem_vars.LWRabsB, LWRabsB{FT,0})
        @test isa(bem_vars.BEMWasteHeat, BEMWasteHeat{FT,0})
        @test isa(bem_vars.BEMEnergyUse, BEMEnergyUse{FT,0})
        @test isa(bem_vars.ParACHeat_ts, ParACHeat_ts{FT,0})
    end

    @testset "TimeSeries" begin
        hours = 24
        tatm = FT(300.0)  # Initial temperature
        qatm = FT(0.01)   # Initial humidity

        bem_vars = initialize_building_energy_model_variables(
            FT, TimeSeries(), hours, tatm, qatm
        )

        # Test structure
        @test bem_vars isa BuildingEnergyModelVariables{FT,1}

        # Test fields initialized to zero
        @test isa(bem_vars.TempVecB, TempVecB{FT,1})
        @test isa(bem_vars.HumidityBuilding, HumidityBuilding{FT,1})
        @test isa(bem_vars.HbuildInt, HbuildInt{FT,1})
        @test isa(bem_vars.LEbuildInt, LEbuildInt{FT,1})
        @test isa(bem_vars.GbuildInt, GbuildInt{FT,1})
        @test isa(bem_vars.SWRabsB, SWRabsB{FT,1})
        @test isa(bem_vars.LWRabsB, LWRabsB{FT,1})
        @test isa(bem_vars.BEMWasteHeat, BEMWasteHeat{FT,1})
        @test isa(bem_vars.BEMEnergyUse, BEMEnergyUse{FT,1})
        @test isa(bem_vars.ParACHeat_ts, ParACHeat_ts{FT,1})

        @test bem_vars.TempVecB.Tceiling[1] == tatm
        @test bem_vars.TempVecB.Tinwallsun[1] == tatm
        @test bem_vars.TempVecB.qbin[1] == qatm
    end
end
