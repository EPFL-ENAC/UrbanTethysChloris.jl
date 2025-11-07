using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Eflux,
    initialize_eflux,
    Runoff,
    initialize_runoff,
    Runon,
    initialize_runon,
    Leakage,
    initialize_leakage,
    Interception,
    initialize_interception,
    dInt_dt,
    initialize_dint_dt,
    Infiltration,
    initialize_infiltration,
    Vwater,
    initialize_vwater,
    dVwater_dt,
    initialize_dvwater_dt,
    Owater,
    initialize_owater,
    OSwater,
    initialize_oswater,
    Qinlat,
    initialize_qinlat,
    SoilPotW,
    initialize_soilpotw,
    ExWater,
    initialize_exwater,
    CiCO2Leaf,
    initialize_cico2leaf,
    WaterFluxVariables,
    initialize_water_flux_variables,
    TimeSlice,
    TimeSeries

using UrbanTethysChloris.ModelComponents.ModelVariables: calculate_soil_values

using UrbanTethysChloris.ModelComponents.Parameters: initialize_parameter_set
using ....TestUtils: load_test_parameters

FT = Float64

input_data = load_test_parameters()

ps = initialize_parameter_set(Float64, input_data)

roof_soil_params = calculate_soil_values(
    ps.soil.roof, ps.vegetation.roof, ps.vegetation.roof
)

ground_soil_params = calculate_soil_values(
    ps.soil.ground, ps.vegetation.tree, ps.vegetation.ground
)

soil_values = (; roof=roof_soil_params, ground=ground_soil_params)

@testset "Water Flux Variables Components" begin
    @testset "Eflux" begin
        @testset "TimeSlice" begin
            eflux = initialize_eflux(FT, TimeSlice())

            # Test structure
            @test eflux isa Eflux{FT,0}

            # Test field access for scalar case
            @test eflux.EfluxRoofImp === 0.0
            @test eflux.EfluxGroundVeg === 0.0
            @test eflux.EfluxUrban === 0.0

            # Test all fields are accessible
            for field in fieldnames(Eflux)
                @test isa(getproperty(eflux, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            eflux = initialize_eflux(FT, TimeSeries(), hours)

            @test eflux isa Eflux{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Eflux)
                @test isa(getproperty(eflux, field), Array{FT,1})
                @test size(getproperty(eflux, field)) == (hours,)
                @test getproperty(eflux, field)[1] == 0
            end
        end
    end

    @testset "Runoff" begin
        @testset "TimeSlice" begin
            runoff = initialize_runoff(FT, TimeSlice())

            # Test structure
            @test runoff isa Runoff{FT,0}

            # Test field access for scalar case
            @test runoff.QRoofImp === 0.0
            @test runoff.QGroundBareSoil === 0.0
            @test runoff.QGroundVegSoil === 0.0

            # Test all fields are accessible
            for field in fieldnames(Runoff)
                @test isa(getproperty(runoff, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            runoff = initialize_runoff(FT, TimeSeries(), hours)

            @test runoff isa Runoff{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Runoff)
                @test isa(getproperty(runoff, field), Array{FT,1})
                @test size(getproperty(runoff, field)) == (hours,)
                @test getproperty(runoff, field)[1] == 0
            end
        end
    end

    @testset "Runon" begin
        @testset "TimeSlice" begin
            runon = initialize_runon(FT, TimeSlice())

            # Test structure
            @test runon isa Runon{FT,0}

            # Test field access for scalar case
            @test runon.RunonRoofTot === 0.0
            @test runon.RunoffGroundTot === 0.0
            @test runon.RunoffUrban === 0.0

            # Test all fields are accessible
            for field in fieldnames(Runon)
                @test isa(getproperty(runon, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            runon = initialize_runon(FT, TimeSeries(), hours)

            @test runon isa Runon{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Runon)
                @test isa(getproperty(runon, field), Array{FT,1})
                @test size(getproperty(runon, field)) == (hours,)
                @test getproperty(runon, field)[1] == 0
            end
        end
    end

    @testset "Leakage" begin
        @testset "TimeSlice" begin
            leakage = initialize_leakage(FT, TimeSlice())

            # Test structure
            @test leakage isa Leakage{FT,0}

            # Test field access for scalar case
            @test leakage.LkRoofImp === 0.0
            @test leakage.LkGroundVeg === 0.0
            @test leakage.LkUrban === 0.0

            # Test all fields are accessible
            for field in fieldnames(Leakage)
                @test isa(getproperty(leakage, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            leakage = initialize_leakage(FT, TimeSeries(), hours)

            @test leakage isa Leakage{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Leakage)
                @test isa(getproperty(leakage, field), Array{FT,1})
                @test size(getproperty(leakage, field)) == (hours,)
                @test getproperty(leakage, field)[1] == 0
            end
        end
    end

    @testset "Interception" begin
        @testset "TimeSlice" begin
            interception = initialize_interception(FT, TimeSlice())

            # Test structure
            @test interception isa Interception{FT,0}

            # Test field access for scalar case
            @test interception.IntRoofImp === 0.0
            @test interception.IntGroundBare === 0.0
            @test interception.IntTree === 0.0

            # Test all fields are accessible
            for field in fieldnames(Interception)
                @test isa(getproperty(interception, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            interception = initialize_interception(FT, TimeSeries(), hours)

            @test interception isa Interception{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Interception)
                @test isa(getproperty(interception, field), Array{FT,1})
                @test size(getproperty(interception, field)) == (hours,)
                @test getproperty(interception, field)[1] == 0
            end
        end
    end

    @testset "dInt_dt" begin
        @testset "TimeSlice" begin
            dint_dt = initialize_dint_dt(FT, TimeSlice())

            # Test structure
            @test dint_dt isa dInt_dt{FT,0}

            # Test field access for scalar case
            @test dint_dt.dInt_dtRoofImp === 0.0
            @test dint_dt.dInt_dtGroundBare === 0.0
            @test dint_dt.dInt_dtTree === 0.0

            # Test all fields are accessible
            for field in fieldnames(dInt_dt)
                @test isa(getproperty(dint_dt, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            dint_dt = initialize_dint_dt(FT, TimeSeries(), hours)

            @test dint_dt isa dInt_dt{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(dInt_dt)
                @test isa(getproperty(dint_dt, field), Array{FT,1})
                @test size(getproperty(dint_dt, field)) == (hours,)
                @test getproperty(dint_dt, field)[1] == 0
            end
        end
    end

    @testset "Infiltration" begin
        @testset "TimeSlice" begin
            infiltration = initialize_infiltration(FT, TimeSlice())

            # Test structure
            @test infiltration isa Infiltration{FT,0}

            # Test field access for scalar case
            @test infiltration.fRoofVeg === 0.0
            @test infiltration.fGroundBare === 0.0
            @test infiltration.fGroundImp === 0.0

            # Test all fields are accessible
            for field in fieldnames(Infiltration)
                @test isa(getproperty(infiltration, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            infiltration = initialize_infiltration(FT, TimeSeries(), hours)

            @test infiltration isa Infiltration{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Infiltration)
                @test isa(getproperty(infiltration, field), Array{FT,1})
                @test size(getproperty(infiltration, field)) == (hours,)
                @test getproperty(infiltration, field)[1] == 0
            end
        end
    end

    @testset "Vwater" begin
        @testset "TimeSlice" begin
            vwater = initialize_vwater(FT, TimeSlice(), soil_values)

            # Test structure
            @test vwater isa Vwater{FT,1}

            # Test field access for scalar case
            @test all(vwater.VRoofSoilVeg .== 0.0)
            @test all(vwater.VGroundSoilBare .== 0.0)
            @test all(vwater.VGroundSoilTot .== 0.0)

            # Test all fields are accessible
            for field in fieldnames(Vwater)
                @test isa(getproperty(vwater, field), Vector{FT})
            end
        end

        @testset "TimeSeries)" begin
            hours = 24
            vwater = initialize_vwater(FT, TimeSeries(), soil_values, hours)

            # Test structure
            @test vwater isa Vwater{FT,2}

            roof_init = soil_values.roof.O33 .* soil_values.roof.dz
            ground_init = soil_values.ground.O33 .* soil_values.ground.dz

            # Test field access for scalar case
            @test vwater.VRoofSoilVeg[1, :] == roof_init
            @test vwater.VGroundSoilBare[1, :] == ground_init
            @test vwater.VGroundSoilTot[1, :] == ground_init

            # Test all fields are accessible
            for field in fieldnames(Vwater)
                @test isa(getproperty(vwater, field), Matrix{FT})
            end
        end
    end

    @testset "dVwater_dt" begin
        @testset "TimeSlice" begin
            dvwater_dt = initialize_dvwater_dt(FT, TimeSlice(), soil_values)

            @test dvwater_dt isa dVwater_dt{FT,1}

            @test size(dvwater_dt.dVRoofSoilVeg_dt) == (soil_values.roof.ms,)
            @test size(dvwater_dt.dVGroundSoilBare_dt) == (soil_values.ground.ms,)
            @test size(dvwater_dt.dVGroundSoilTot_dt) == (soil_values.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(dVwater_dt)
                @test isa(getproperty(dvwater_dt, field), Array{FT,1})
                @test getproperty(dvwater_dt, field)[1] == 0
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            dvwater_dt = initialize_dvwater_dt(FT, TimeSeries(), soil_values, hours)

            # Test structure
            @test dvwater_dt isa dVwater_dt{FT,2}

            @test size(dvwater_dt.dVRoofSoilVeg_dt) == (hours, soil_values.roof.ms)
            @test size(dvwater_dt.dVGroundSoilBare_dt) == (hours, soil_values.ground.ms)
            @test size(dvwater_dt.dVGroundSoilTot_dt) == (hours, soil_values.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(dVwater_dt)
                @test isa(getproperty(dvwater_dt, field), Matrix{FT})
            end
        end
    end

    @testset "Owater" begin
        @testset "TimeSlice" begin
            owater = initialize_owater(FT, TimeSlice(), soil_values)

            @test owater isa Owater{FT,1}

            @test size(owater.OwRoofSoilVeg) == (soil_values.roof.ms,)
            @test size(owater.OwGroundSoilBare) == (soil_values.ground.ms,)
            @test size(owater.OwGroundSoilTot) == (soil_values.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Owater)
                @test isa(getproperty(owater, field), Array{FT,1})
                @test getproperty(owater, field)[1] == 0
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            owater = initialize_owater(FT, TimeSeries(), soil_values, hours)

            # Test structure
            @test owater isa Owater{FT,2}

            @test size(owater.OwRoofSoilVeg) == (hours, soil_values.roof.ms)
            @test size(owater.OwGroundSoilBare) == (hours, soil_values.ground.ms)
            @test size(owater.OwGroundSoilTot) == (hours, soil_values.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(Owater)
                @test isa(getproperty(owater, field), Matrix{FT})
            end
        end
    end
    @testset "OSwater" begin
        @testset "TimeSlice" begin
            oswater = initialize_oswater(FT, TimeSlice(), soil_values)

            @test oswater isa OSwater{FT,1}

            @test size(oswater.OSwRoofSoilVeg) == (soil_values.roof.ms,)
            @test size(oswater.OSwGroundSoilBare) == (soil_values.ground.ms,)
            @test size(oswater.OSwGroundSoilTot) == (soil_values.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(OSwater)
                @test isa(getproperty(oswater, field), Array{FT,1})
                @test getproperty(oswater, field)[1] == 0
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            oswater = initialize_oswater(FT, TimeSeries(), soil_values, hours)

            # Test structure
            @test oswater isa OSwater{FT,2}

            @test size(oswater.OSwRoofSoilVeg) == (hours, soil_values.roof.ms)
            @test size(oswater.OSwGroundSoilBare) == (hours, soil_values.ground.ms)
            @test size(oswater.OSwGroundSoilTot) == (hours, soil_values.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(OSwater)
                @test isa(getproperty(oswater, field), Matrix{FT})
            end
        end
    end

    @testset "Qinlat" begin
        @testset "TimeSlice" begin
            qinlat = initialize_qinlat(FT, TimeSlice(), soil_values)

            @test qinlat isa Qinlat{FT,1}

            @test size(qinlat.Qin_bare2imp) == (soil_values.ground.ms,)
            @test size(qinlat.Qin_imp2bare) == (soil_values.ground.ms,)
            @test size(qinlat.Qin_imp) == (soil_values.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Qinlat)
                @test isa(getproperty(qinlat, field), Array{FT,1})
                @test getproperty(qinlat, field)[1] == 0
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            qinlat = initialize_qinlat(FT, TimeSeries(), soil_values, hours)

            # Test structure
            @test qinlat isa Qinlat{FT,2}

            @test size(qinlat.Qin_bare2imp) == (hours, soil_values.ground.ms)
            @test size(qinlat.Qin_imp2bare) == (hours, soil_values.ground.ms)
            @test size(qinlat.Qin_imp) == (hours, soil_values.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(Qinlat)
                @test isa(getproperty(qinlat, field), Matrix{FT})
            end
        end
    end

    @testset "ExWater" begin
        @testset "TimeSlice" begin
            exwater = initialize_exwater(FT, TimeSlice(), soil_values)

            @test exwater isa ExWater{FT,1}

            @test size(exwater.ExWaterRoofVeg_H) == (soil_values.roof.ms,)
            @test size(exwater.ExWaterRoofVeg_L) == (soil_values.roof.ms,)
            @test size(exwater.ExWaterGroundImp_H) == (soil_values.ground.ms,)
            @test size(exwater.ExWaterGroundImp_L) == (soil_values.ground.ms,)
            @test size(exwater.ExWaterGroundBare_H) == (soil_values.ground.ms,)
            @test size(exwater.ExWaterGroundBare_L) == (soil_values.ground.ms,)
            @test size(exwater.ExWaterGroundVeg_H) == (soil_values.ground.ms,)
            @test size(exwater.ExWaterGroundVeg_L) == (soil_values.ground.ms,)
            @test size(exwater.ExWaterGroundTot_H) == (soil_values.ground.ms,)
            @test size(exwater.ExWaterGroundTot_L) == (soil_values.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(ExWater)
                @test isa(getproperty(exwater, field), Array{FT,1})
                @test getproperty(exwater, field)[1] == 0
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            exwater = initialize_exwater(FT, TimeSeries(), soil_values, hours)

            # Test structure
            @test exwater isa ExWater{FT,2}

            @test size(exwater.ExWaterRoofVeg_H) == (hours, soil_values.roof.ms)
            @test size(exwater.ExWaterRoofVeg_L) == (hours, soil_values.roof.ms)
            @test size(exwater.ExWaterGroundImp_H) == (hours, soil_values.ground.ms)
            @test size(exwater.ExWaterGroundImp_L) == (hours, soil_values.ground.ms)
            @test size(exwater.ExWaterGroundBare_H) == (hours, soil_values.ground.ms)
            @test size(exwater.ExWaterGroundBare_L) == (hours, soil_values.ground.ms)
            @test size(exwater.ExWaterGroundVeg_H) == (hours, soil_values.ground.ms)
            @test size(exwater.ExWaterGroundVeg_L) == (hours, soil_values.ground.ms)
            @test size(exwater.ExWaterGroundTot_H) == (hours, soil_values.ground.ms)
            @test size(exwater.ExWaterGroundTot_L) == (hours, soil_values.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(ExWater)
                @test isa(getproperty(exwater, field), Matrix{FT})
            end
        end
    end

    @testset "SoilPotW" begin
        @testset "TimeSlice" begin
            soilpotw = initialize_soilpotw(FT, TimeSlice())

            # Test structure
            @test soilpotw isa SoilPotW{FT,0}

            # Test field access for scalar case
            @test soilpotw.SoilPotWRoofVeg_H === 0.0
            @test soilpotw.SoilPotWGroundImp_H === 0.0
            @test soilpotw.SoilPotWGroundVeg_L === 0.0

            # Test all fields are accessible
            for field in fieldnames(SoilPotW)
                @test isa(getproperty(soilpotw, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            soilpotw = initialize_soilpotw(FT, TimeSeries(), hours)

            @test soilpotw isa SoilPotW{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(SoilPotW)
                @test isa(getproperty(soilpotw, field), Array{FT,1})
                @test size(getproperty(soilpotw, field)) == (hours,)
                @test getproperty(soilpotw, field)[1] == 0
            end
        end
    end

    @testset "CiCO2Leaf" begin
        @testset "TimeSlice" begin
            cico2leaf = initialize_cico2leaf(FT, TimeSlice())

            # Test structure
            @test cico2leaf isa CiCO2Leaf{FT,0}

            # Test field access for scalar case
            @test cico2leaf.CiCO2LeafRoofVegSun === 0.0
            @test cico2leaf.CiCO2LeafGroundVegSun === 0.0
            @test cico2leaf.CiCO2LeafTreeSun === 0.0

            # Test all fields are accessible
            for field in fieldnames(CiCO2Leaf)
                @test isa(getproperty(cico2leaf, field), FT)
            end
        end

        @testset "TimeSeries" begin
            hours = 24
            initial_value = 300.0
            cico2leaf = initialize_cico2leaf(FT, TimeSeries(), initial_value, hours)

            @test cico2leaf isa CiCO2Leaf{FT,1}

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(CiCO2Leaf)
                @test isa(getproperty(cico2leaf, field), Array{FT,1})
                @test size(getproperty(cico2leaf, field)) == (hours,)
                @test getproperty(cico2leaf, field)[1] == initial_value
            end
        end
    end
end

@testset "WaterFluxVariables" begin
    @testset "TimeSlice" begin
        water_flux_vars = initialize_water_flux_variables(
            FT, TimeSlice(), ps.soil, ps.vegetation
        )

        # Test structure
        @test water_flux_vars isa WaterFluxVariables{FT,0,1}

        # Test fields are properly initialized
        @test isa(water_flux_vars.Eflux, Eflux{FT,0})
        @test isa(water_flux_vars.Runoff, Runoff{FT,0})
        @test isa(water_flux_vars.Runon, Runon{FT,0})
        @test isa(water_flux_vars.Leakage, Leakage{FT,0})
        @test isa(water_flux_vars.Interception, Interception{FT,0})
        @test isa(water_flux_vars.dInt_dt, dInt_dt{FT,0})
        @test isa(water_flux_vars.Infiltration, Infiltration{FT,0})
        @test isa(water_flux_vars.SoilPotW, SoilPotW{FT,0})

        # Soil layer fields should be properly initialized with correct dimensions
        @test isa(water_flux_vars.Vwater, Vwater{FT,1})
        @test isa(water_flux_vars.dVwater_dt, dVwater_dt{FT,1})
        @test isa(water_flux_vars.Owater, Owater{FT,1})
        @test isa(water_flux_vars.OSwater, OSwater{FT,1})
        @test isa(water_flux_vars.Qinlat, Qinlat{FT,1})
        @test isa(water_flux_vars.ExWater, ExWater{FT,1})
    end

    @testset "TimeSeries" begin
        # Mock the water_flux_variables initialization
        hours = 24
        initial_value = 300.0
        water_flux_vars = initialize_water_flux_variables(
            FT, TimeSeries(), ps.soil, ps.vegetation, initial_value, hours
        )

        # Test structure
        @test water_flux_vars isa WaterFluxVariables{FT,1,2}

        # Test fields are properly initialized
        @test isa(water_flux_vars.Eflux, Eflux{FT,1})
        @test isa(water_flux_vars.Runoff, Runoff{FT,1})
        @test isa(water_flux_vars.Runon, Runon{FT,1})
        @test isa(water_flux_vars.Leakage, Leakage{FT,1})
        @test isa(water_flux_vars.Interception, Interception{FT,1})
        @test isa(water_flux_vars.dInt_dt, dInt_dt{FT,1})
        @test isa(water_flux_vars.Infiltration, Infiltration{FT,1})
        @test isa(water_flux_vars.SoilPotW, SoilPotW{FT,1})

        # Soil layer fields should be properly initialized with correct dimensions
        @test isa(water_flux_vars.Vwater, Vwater{FT,2})
        @test isa(water_flux_vars.dVwater_dt, dVwater_dt{FT,2})
        @test isa(water_flux_vars.Owater, Owater{FT,2})
        @test isa(water_flux_vars.OSwater, OSwater{FT,2})
        @test isa(water_flux_vars.Qinlat, Qinlat{FT,2})
        @test isa(water_flux_vars.ExWater, ExWater{FT,2})
    end
end

@testset "get/setindex" begin
    hours = 24
    initial_value = 300.0
    water_flux_vars = initialize_water_flux_variables(
        FT, TimeSeries(), ps.soil, ps.vegetation, initial_value, hours
    )

    EfluxRoofImp = FT(1.0)
    QRoofImp = FT(2.0)
    RunonRoofTot = FT(3.0)
    LkRoofImp = FT(4.0)
    IntRoofImp = FT(5.0)
    dInt_dtRoofImp = FT(6.0)
    fRoofVeg = FT(7.0)
    VRoofSoilVeg = fill(FT(8.0), ps.soil.roof.ms)
    dVGroundSoilVeg_dt = fill(FT(9.0), ps.soil.ground.ms)
    OwRoofSoilVeg = fill(FT(10.0), ps.soil.roof.ms)
    OSwRoofSoilVeg = fill(FT(11.0), ps.soil.roof.ms)
    Qin_bare2imp = fill(FT(12.0), ps.soil.ground.ms)
    ExWaterRoofVeg_H = fill(FT(13.0), ps.soil.roof.ms)
    SoilPotWRoofVeg_H = FT(14.0)
    CiCO2LeafRoofVegSun = FT(15.0)

    x = water_flux_vars[1]
    x.Eflux.EfluxRoofImp = EfluxRoofImp
    x.Runoff.QRoofImp = QRoofImp
    x.Runon.RunonRoofTot = RunonRoofTot
    x.Leakage.LkRoofImp = LkRoofImp
    x.Interception.IntRoofImp = IntRoofImp
    x.dInt_dt.dInt_dtRoofImp = dInt_dtRoofImp
    x.Infiltration.fRoofVeg = fRoofVeg
    x.Vwater.VRoofSoilVeg = VRoofSoilVeg
    x.dVwater_dt.dVGroundSoilVeg_dt = dVGroundSoilVeg_dt
    x.Owater.OwRoofSoilVeg = OwRoofSoilVeg
    x.OSwater.OSwRoofSoilVeg = OSwRoofSoilVeg
    x.Qinlat.Qin_bare2imp = Qin_bare2imp
    x.ExWater.ExWaterRoofVeg_H = ExWaterRoofVeg_H
    x.SoilPotW.SoilPotWRoofVeg_H = SoilPotWRoofVeg_H
    x.CiCO2Leaf.CiCO2LeafRoofVegSun = CiCO2LeafRoofVegSun

    water_flux_vars[2] = x
    @test water_flux_vars.Eflux.EfluxRoofImp[2] == EfluxRoofImp
    @test water_flux_vars.Runoff.QRoofImp[2] == QRoofImp
    @test water_flux_vars.Runon.RunonRoofTot[2] == RunonRoofTot
    @test water_flux_vars.Leakage.LkRoofImp[2] == LkRoofImp
    @test water_flux_vars.Interception.IntRoofImp[2] == IntRoofImp
    @test water_flux_vars.dInt_dt.dInt_dtRoofImp[2] == dInt_dtRoofImp
    @test water_flux_vars.Infiltration.fRoofVeg[2] == fRoofVeg
    @test water_flux_vars.Vwater.VRoofSoilVeg[2, :] == VRoofSoilVeg
    @test water_flux_vars.dVwater_dt.dVGroundSoilVeg_dt[2, :] == dVGroundSoilVeg_dt
    @test water_flux_vars.Owater.OwRoofSoilVeg[2, :] == OwRoofSoilVeg
    @test water_flux_vars.OSwater.OSwRoofSoilVeg[2, :] == OSwRoofSoilVeg
    @test water_flux_vars.Qinlat.Qin_bare2imp[2, :] == Qin_bare2imp
    @test water_flux_vars.ExWater.ExWaterRoofVeg_H[2, :] == ExWaterRoofVeg_H
    @test water_flux_vars.SoilPotW.SoilPotWRoofVeg_H[2] == SoilPotWRoofVeg_H
    @test water_flux_vars.CiCO2Leaf.CiCO2LeafRoofVegSun[2] == CiCO2LeafRoofVegSun
end
