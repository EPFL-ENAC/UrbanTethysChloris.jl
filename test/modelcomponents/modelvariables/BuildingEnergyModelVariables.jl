using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    TempVecB,
    HumidityBuilding,
    HbuildInt,
    LEbuildInt,
    GbuildInt,
    SWRabsB,
    LWRabsB,
    BEMWasteHeat,
    BEMEnergyUse,
    ParACHeat_ts,
    BuildingEnergyModelVariables

FT = Float64

@testset "Subsets" begin
    @testset "TempVecB" begin
        tempvecb = TempVecB(FT)

        # Test structure
        @test tempvecb isa TempVecB{FT}

        # Test all fields are accessible
        for field in fieldnames(TempVecB)
            @test isa(getproperty(tempvecb, field), FT)
            @test getproperty(tempvecb, field) == 0
        end
    end

    @testset "HumidityBuilding" begin
        humidityb = HumidityBuilding(FT)

        # Test structure
        @test humidityb isa HumidityBuilding{FT}

        # Test all fields are accessible
        for field in fieldnames(HumidityBuilding)
            @test isa(getproperty(humidityb, field), FT)
            @test getproperty(humidityb, field) == 0
        end
    end

    @testset "HbuildInt" begin
        hbuildint = HbuildInt(FT)

        # Test structure
        @test hbuildint isa HbuildInt{FT}

        # Test all fields are accessible
        for field in fieldnames(HbuildInt)
            @test isa(getproperty(hbuildint, field), FT)
            @test getproperty(hbuildint, field) == 0
        end
    end

    @testset "LEbuildInt" begin
        lebuildint = LEbuildInt(FT)

        # Test structure
        @test lebuildint isa LEbuildInt{FT}

        # Test all fields are accessible
        for field in fieldnames(LEbuildInt)
            @test isa(getproperty(lebuildint, field), FT)
            @test getproperty(lebuildint, field) == 0
        end
    end

    @testset "GbuildInt" begin
        gbuildint = GbuildInt(FT)

        # Test structure
        @test gbuildint isa GbuildInt{FT}

        # Test all fields are accessible
        for field in fieldnames(GbuildInt)
            @test isa(getproperty(gbuildint, field), FT)
            @test getproperty(gbuildint, field) == 0
        end
    end

    @testset "SWRabsB" begin
        swrabsb = SWRabsB(FT)

        # Test structure
        @test swrabsb isa SWRabsB{FT}

        # Test all fields are accessible
        for field in fieldnames(SWRabsB)
            @test isa(getproperty(swrabsb, field), FT)
            @test getproperty(swrabsb, field) == 0
        end
    end

    @testset "LWRabsB" begin
        lwrabsb = LWRabsB(FT)

        # Test structure
        @test lwrabsb isa LWRabsB{FT}

        # Test all fields are accessible
        for field in fieldnames(LWRabsB)
            @test isa(getproperty(lwrabsb, field), FT)
            @test getproperty(lwrabsb, field) == 0
        end
    end

    @testset "BEMWasteHeat" begin
        bemwasteheat = BEMWasteHeat(FT)

        # Test structure
        @test bemwasteheat isa BEMWasteHeat{FT}

        # Test all fields are accessible
        for field in fieldnames(BEMWasteHeat)
            @test isa(getproperty(bemwasteheat, field), FT)
            @test getproperty(bemwasteheat, field) == 0
        end
    end

    @testset "BEMEnergyUse" begin
        bemenergyuse = BEMEnergyUse(FT)

        # Test structure
        @test bemenergyuse isa BEMEnergyUse{FT}

        # Test all fields are accessible
        for field in fieldnames(BEMEnergyUse)
            @test isa(getproperty(bemenergyuse, field), FT)
            @test getproperty(bemenergyuse, field) == 0
        end
    end

    @testset "ParACHeat_ts" begin
        paracheat_ts = ParACHeat_ts(FT)

        # Test structure
        @test paracheat_ts isa ParACHeat_ts{FT}

        # Test all fields are accessible
        for field in fieldnames(ParACHeat_ts)
            @test isa(getproperty(paracheat_ts, field), FT)
            @test getproperty(paracheat_ts, field) == 0
        end
    end
end

@testset "BuildingEnergyModelVariables" begin
    bem_vars = BuildingEnergyModelVariables(FT)

    # Test structure
    @test bem_vars isa BuildingEnergyModelVariables{FT}

    # Test fields initialized to zero
    @test isa(bem_vars.TempVecB, TempVecB{FT})
    @test isa(bem_vars.HumidityBuilding, HumidityBuilding{FT})
    @test isa(bem_vars.HbuildInt, HbuildInt{FT})
    @test isa(bem_vars.LEbuildInt, LEbuildInt{FT})
    @test isa(bem_vars.GbuildInt, GbuildInt{FT})
    @test isa(bem_vars.SWRabsB, SWRabsB{FT})
    @test isa(bem_vars.LWRabsB, LWRabsB{FT})
    @test isa(bem_vars.BEMWasteHeat, BEMWasteHeat{FT})
    @test isa(bem_vars.BEMEnergyUse, BEMEnergyUse{FT})
    @test isa(bem_vars.ParACHeat_ts, ParACHeat_ts{FT})
end
