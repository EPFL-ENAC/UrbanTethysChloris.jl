using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Eflux,
    Runoff,
    Runon,
    Leakage,
    Interception,
    dInt_dt,
    Infiltration,
    Vwater,
    dVwater_dt,
    Owater,
    OSwater,
    Qinlat,
    SoilPotW,
    ExWater,
    CiCO2Leaf,
    WaterFluxVariables,
    ground_fields,
    roof_fields
using StaticArrays

using UrbanTethysChloris.ModelComponents.Parameters: initialize_soil_parameters
using ....TestUtils: load_test_parameters

FT = Float64

input_data = load_test_parameters()

soil = initialize_soil_parameters(Float64, input_data["soil"])

@testset "Water Flux Variables Components" begin
    @testset "Eflux" begin
        eflux = Eflux(FT)

        # Test structure
        @test eflux isa Eflux{FT}

        # Test all fields are accessible
        for field in fieldnames(Eflux)
            @test isa(getproperty(eflux, field), FT)
            @test getproperty(eflux, field) == 0
        end
    end

    @testset "Runoff" begin
        runoff = Runoff(FT)

        # Test structure
        @test runoff isa Runoff{FT}

        # Test all fields are accessible
        for field in fieldnames(Runoff)
            @test isa(getproperty(runoff, field), FT)
            @test getproperty(runoff, field) == 0
        end
    end

    @testset "Runon" begin
        runon = Runon(FT)

        # Test structure
        @test runon isa Runon{FT}

        # Test all fields are accessible
        for field in fieldnames(Runon)
            @test isa(getproperty(runon, field), FT)
            @test getproperty(runon, field) == 0
        end
    end

    @testset "Leakage" begin
        leakage = Leakage(FT)

        # Test structure
        @test leakage isa Leakage{FT}

        # Test all fields are accessible
        for field in fieldnames(Leakage)
            @test isa(getproperty(leakage, field), FT)
            @test getproperty(leakage, field) == 0
        end
    end

    @testset "Interception" begin
        interception = Interception(FT)

        # Test structure
        @test interception isa Interception{FT}

        # Test all fields are accessible
        for field in fieldnames(Interception)
            @test isa(getproperty(interception, field), FT)
            @test getproperty(interception, field) == 0
        end
    end

    @testset "dInt_dt" begin
        dint_dt = dInt_dt(FT)

        # Test structure
        @test dint_dt isa dInt_dt{FT}

        # Test all fields are accessible
        for field in fieldnames(dInt_dt)
            @test isa(getproperty(dint_dt, field), FT)
            @test getproperty(dint_dt, field) == 0
        end
    end

    @testset "Infiltration" begin
        infiltration = Infiltration(FT)

        # Test structure
        @test infiltration isa Infiltration{FT}

        # Test all fields are accessible
        for field in fieldnames(Infiltration)
            @test isa(getproperty(infiltration, field), FT)
            @test getproperty(infiltration, field) == 0
        end
    end

    @testset "Vwater" begin
        vwater = Vwater(FT, soil)

        # Test structure
        @test vwater isa Vwater{FT,soil.roof.ms,soil.ground.ms}

        # Test all fields are accessible
        for field in Symbol.(ground_fields(Vwater))
            @test isa(getproperty(vwater, field), Vector{FT})
            @test all(getproperty(vwater, field) .== 0)
        end

        for field in Symbol.(roof_fields(Vwater))
            @test isa(getproperty(vwater, field), Vector{FT})
            @test all(getproperty(vwater, field) .== 0)
        end
    end

    @testset "dVwater_dt" begin
        dvwater_dt = dVwater_dt(FT)

        @test dvwater_dt isa dVwater_dt{FT}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in Symbol.(ground_fields(dVwater_dt))
            @test isa(getproperty(dvwater_dt, field), FT)
            @test all(getproperty(dvwater_dt, field) == 0)
        end

        for field in Symbol.(roof_fields(dVwater_dt))
            @test isa(getproperty(dvwater_dt, field), FT)
            @test all(getproperty(dvwater_dt, field) == 0)
        end
    end

    @testset "Owater" begin
        owater = Owater(FT, soil)

        @test owater isa Owater{FT,soil.roof.ms,soil.ground.ms}

        # Test all fields are accessible
        for field in Symbol.(ground_fields(Owater))
            @test isa(getproperty(owater, field), Vector{FT})
            @test all(getproperty(owater, field) .== 0)
        end

        for field in Symbol.(roof_fields(Owater))
            @test isa(getproperty(owater, field), Vector{FT})
            @test all(getproperty(owater, field) .== 0)
        end
    end
    @testset "OSwater" begin
        oswater = OSwater(FT)

        @test oswater isa OSwater{FT}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in Symbol.(ground_fields(OSwater))
            @test isa(getproperty(oswater, field), FT)
            @test all(getproperty(oswater, field) == 0)
        end

        for field in Symbol.(roof_fields(OSwater))
            @test isa(getproperty(oswater, field), FT)
            @test all(getproperty(oswater, field) == 0)
        end
    end

    @testset "Qinlat" begin
        qinlat = Qinlat(FT, soil)

        @test qinlat isa Qinlat{FT,soil.ground.ms}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in Symbol.(ground_fields(Qinlat))
            @test isa(getproperty(qinlat, field), Vector{FT})
            @test all(getproperty(qinlat, field) .== 0)
        end
    end

    @testset "ExWater" begin
        exwater = ExWater(FT, soil)

        @test exwater isa ExWater{FT,soil.roof.ms,soil.ground.ms}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in Symbol.(ground_fields(ExWater))
            @test isa(getproperty(exwater, field), Vector{FT})
            @test all(getproperty(exwater, field) .== 0)
        end

        for field in Symbol.(roof_fields(ExWater))
            @test isa(getproperty(exwater, field), Vector{FT})
            @test all(getproperty(exwater, field) .== 0)
        end
    end

    @testset "SoilPotW" begin
        soilpotw = SoilPotW(FT)

        # Test structure
        @test soilpotw isa SoilPotW{FT}

        # Test all fields are accessible
        for field in fieldnames(SoilPotW)
            @test isa(getproperty(soilpotw, field), FT)
            @test getproperty(soilpotw, field) == 0
        end
    end

    @testset "CiCO2Leaf" begin
        cico2leaf = CiCO2Leaf(FT)

        # Test structure
        @test cico2leaf isa CiCO2Leaf{FT}

        # Test all fields are accessible
        for field in fieldnames(CiCO2Leaf)
            @test isa(getproperty(cico2leaf, field), FT)
            @test getproperty(cico2leaf, field) == 0
        end
    end
end

@testset "WaterFluxVariables" begin
    MR = soil.roof.ms
    MG = soil.ground.ms
    water_flux_vars = WaterFluxVariables(FT, soil)

    # Test structure
    @test water_flux_vars isa WaterFluxVariables{FT,MR,MG}

    # Test fields are properly initialized
    @test isa(water_flux_vars.Eflux, Eflux{FT})
    @test isa(water_flux_vars.Runoff, Runoff{FT})
    @test isa(water_flux_vars.Runon, Runon{FT})
    @test isa(water_flux_vars.Leakage, Leakage{FT})
    @test isa(water_flux_vars.Interception, Interception{FT})
    @test isa(water_flux_vars.dInt_dt, dInt_dt{FT})
    @test isa(water_flux_vars.Infiltration, Infiltration{FT})
    @test isa(water_flux_vars.SoilPotW, SoilPotW{FT})

    # Soil layer fields should be properly initialized with correct dimensions
    @test isa(water_flux_vars.Vwater, Vwater{FT,MR,MG})
    @test isa(water_flux_vars.dVwater_dt, dVwater_dt{FT})
    @test isa(water_flux_vars.Owater, Owater{FT,MR,MG})
    @test isa(water_flux_vars.OSwater, OSwater{FT})
    @test isa(water_flux_vars.Qinlat, Qinlat{FT,MG})
    @test isa(water_flux_vars.ExWater, ExWater{FT,MR,MG})
end
