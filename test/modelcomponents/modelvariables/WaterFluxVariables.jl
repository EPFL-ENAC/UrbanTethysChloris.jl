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
    initialize_water_flux_variables

using UrbanTethysChloris.ModelComponents.Parameters: initialize_parameter_set
using ....TestUtils: load_test_parameters

FT = Float64

input_data = load_test_parameters()

ps = initialize_parameter_set(Float64, input_data)

@testset "Water Flux Variables Components" begin
    @testset "Eflux scalar (N=0)" begin
        eflux = initialize_eflux(FT, 0)

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

    @testset "Eflux vector (N=1)" begin
        hours = 24
        eflux = initialize_eflux(FT, 1, hours)

        @test eflux isa Eflux{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Eflux)
            @test isa(getproperty(eflux, field), Array{FT,1})
            @test size(getproperty(eflux, field)) == (hours,)
            @test getproperty(eflux, field)[1] == 0
        end
    end

    @testset "Runoff scalar (N=0)" begin
        runoff = initialize_runoff(FT, 0)

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

    @testset "Runoff vector (N=1)" begin
        hours = 24
        runoff = initialize_runoff(FT, 1, hours)

        @test runoff isa Runoff{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Runoff)
            @test isa(getproperty(runoff, field), Array{FT,1})
            @test size(getproperty(runoff, field)) == (hours,)
            @test getproperty(runoff, field)[1] == 0
        end
    end

    @testset "Runon scalar (N=0)" begin
        runon = initialize_runon(FT, 0)

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

    @testset "Runon vector (N=1)" begin
        hours = 24
        runon = initialize_runon(FT, 1, hours)

        @test runon isa Runon{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Runon)
            @test isa(getproperty(runon, field), Array{FT,1})
            @test size(getproperty(runon, field)) == (hours,)
            @test getproperty(runon, field)[1] == 0
        end
    end

    @testset "Leakage scalar (N=0)" begin
        leakage = initialize_leakage(FT, 0)

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

    @testset "Leakage vector (N=1)" begin
        hours = 24
        leakage = initialize_leakage(FT, 1, hours)

        @test leakage isa Leakage{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Leakage)
            @test isa(getproperty(leakage, field), Array{FT,1})
            @test size(getproperty(leakage, field)) == (hours,)
            @test getproperty(leakage, field)[1] == 0
        end
    end

    @testset "Interception scalar (N=0)" begin
        interception = initialize_interception(FT, 0)

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

    @testset "Interception vector (N=1)" begin
        hours = 24
        interception = initialize_interception(FT, 1, hours)

        @test interception isa Interception{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Interception)
            @test isa(getproperty(interception, field), Array{FT,1})
            @test size(getproperty(interception, field)) == (hours,)
            @test getproperty(interception, field)[1] == 0
        end
    end

    @testset "dInt_dt scalar (N=0)" begin
        dint_dt = initialize_dint_dt(FT, 0)

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

    @testset "dInt_dt vector (N=1)" begin
        hours = 24
        dint_dt = initialize_dint_dt(FT, 1, hours)

        @test dint_dt isa dInt_dt{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(dInt_dt)
            @test isa(getproperty(dint_dt, field), Array{FT,1})
            @test size(getproperty(dint_dt, field)) == (hours,)
            @test getproperty(dint_dt, field)[1] == 0
        end
    end

    @testset "Infiltration scalar (N=0)" begin
        infiltration = initialize_infiltration(FT, 0)

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

    @testset "Infiltration vector (N=1)" begin
        hours = 24
        infiltration = initialize_infiltration(FT, 1, hours)

        @test infiltration isa Infiltration{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Infiltration)
            @test isa(getproperty(infiltration, field), Array{FT,1})
            @test size(getproperty(infiltration, field)) == (hours,)
            @test getproperty(infiltration, field)[1] == 0
        end
    end

    @testset "Vwater vector (N=1)" begin
        vwater = initialize_vwater(FT, 1, ps.soil)

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

    @testset "Vwater matrix (N=2)" begin
        hours = 24
        vwater = initialize_vwater(FT, 2, ps.soil, hours)

        # Test structure
        @test vwater isa Vwater{FT,2}

        roof_init = ps.soil.roof.O33 .* ps.soil.roof.dz
        ground_init = ps.soil.ground.O33 .* ps.soil.ground.dz

        # Test field access for scalar case
        @test vwater.VRoofSoilVeg[1, :] == roof_init
        @test vwater.VGroundSoilBare[1, :] == ground_init
        @test vwater.VGroundSoilTot[1, :] == ground_init

        # Test all fields are accessible
        for field in fieldnames(Vwater)
            @test isa(getproperty(vwater, field), Matrix{FT})
        end
    end

    @testset "dVwater_dt vector (N=1)" begin
        dvwater_dt = initialize_dvwater_dt(FT, 1, ps.soil)

        @test dvwater_dt isa dVwater_dt{FT,1}

        @test size(dvwater_dt.dVRoofSoilVeg_dt) == (ps.soil.roof.ms,)
        @test size(dvwater_dt.dVGroundSoilBare_dt) == (ps.soil.ground.ms,)
        @test size(dvwater_dt.dVGroundSoilTot_dt) == (ps.soil.ground.ms,)

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(dVwater_dt)
            @test isa(getproperty(dvwater_dt, field), Array{FT,1})
            @test getproperty(dvwater_dt, field)[1] == 0
        end
    end

    @testset "dVwater_dt matrix (N=2)" begin
        hours = 24
        dvwater_dt = initialize_dvwater_dt(FT, 2, ps.soil, hours)

        # Test structure
        @test dvwater_dt isa dVwater_dt{FT,2}

        @test size(dvwater_dt.dVRoofSoilVeg_dt) == (hours, ps.soil.roof.ms)
        @test size(dvwater_dt.dVGroundSoilBare_dt) == (hours, ps.soil.ground.ms)
        @test size(dvwater_dt.dVGroundSoilTot_dt) == (hours, ps.soil.ground.ms)

        # Test all fields are accessible
        for field in fieldnames(dVwater_dt)
            @test isa(getproperty(dvwater_dt, field), Matrix{FT})
        end
    end

    @testset "Owater vector (N=1)" begin
        owater = initialize_owater(FT, 1, ps.soil)

        @test owater isa Owater{FT,1}

        @test size(owater.OwRoofSoilVeg) == (ps.soil.roof.ms,)
        @test size(owater.OwGroundSoilBare) == (ps.soil.ground.ms,)
        @test size(owater.OwGroundSoilTot) == (ps.soil.ground.ms,)

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Owater)
            @test isa(getproperty(owater, field), Array{FT,1})
            @test getproperty(owater, field)[1] == 0
        end
    end

    @testset "Owater matrix (N=2)" begin
        hours = 24
        owater = initialize_owater(FT, 2, ps.soil, hours)

        # Test structure
        @test owater isa Owater{FT,2}

        @test size(owater.OwRoofSoilVeg) == (hours, ps.soil.roof.ms)
        @test size(owater.OwGroundSoilBare) == (hours, ps.soil.ground.ms)
        @test size(owater.OwGroundSoilTot) == (hours, ps.soil.ground.ms)

        # Test all fields are accessible
        for field in fieldnames(Owater)
            @test isa(getproperty(owater, field), Matrix{FT})
        end
    end

    @testset "OSwater" begin
        @testset "OSwater vector (N=1)" begin
            oswater = initialize_oswater(FT, 1, ps.soil)

            @test oswater isa OSwater{FT,1}

            @test size(oswater.OSwRoofSoilVeg) == (ps.soil.roof.ms,)
            @test size(oswater.OSwGroundSoilBare) == (ps.soil.ground.ms,)
            @test size(oswater.OSwGroundSoilTot) == (ps.soil.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(OSwater)
                @test isa(getproperty(oswater, field), Array{FT,1})
                @test getproperty(oswater, field)[1] == 0
            end
        end

        @testset "OSwater matrix (N=2)" begin
            hours = 24
            oswater = initialize_oswater(FT, 2, ps.soil, hours)

            # Test structure
            @test oswater isa OSwater{FT,2}

            @test size(oswater.OSwRoofSoilVeg) == (hours, ps.soil.roof.ms)
            @test size(oswater.OSwGroundSoilBare) == (hours, ps.soil.ground.ms)
            @test size(oswater.OSwGroundSoilTot) == (hours, ps.soil.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(OSwater)
                @test isa(getproperty(oswater, field), Matrix{FT})
            end
        end
    end

    @testset "Qinlat" begin
        @testset "Qinlat vector (N=1)" begin
            qinlat = initialize_qinlat(FT, 1, ps.soil)

            @test qinlat isa Qinlat{FT,1}

            @test size(qinlat.Qin_bare2imp) == (ps.soil.ground.ms,)
            @test size(qinlat.Qin_imp2bare) == (ps.soil.ground.ms,)
            @test size(qinlat.Qin_imp) == (ps.soil.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(Qinlat)
                @test isa(getproperty(qinlat, field), Array{FT,1})
                @test getproperty(qinlat, field)[1] == 0
            end
        end

        @testset "Qinlat matrix (N=2)" begin
            hours = 24
            qinlat = initialize_qinlat(FT, 2, ps.soil, hours)

            # Test structure
            @test qinlat isa Qinlat{FT,2}

            @test size(qinlat.Qin_bare2imp) == (hours, ps.soil.ground.ms)
            @test size(qinlat.Qin_imp2bare) == (hours, ps.soil.ground.ms)
            @test size(qinlat.Qin_imp) == (hours, ps.soil.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(Qinlat)
                @test isa(getproperty(qinlat, field), Matrix{FT})
            end
        end
    end

    @testset "ExWater" begin
        @testset "ExWater vector (N=1)" begin
            exwater = initialize_exwater(FT, 1, ps.soil)

            @test exwater isa ExWater{FT,1}

            @test size(exwater.ExWaterRoofVeg_H) == (ps.soil.roof.ms,)
            @test size(exwater.ExWaterRoofVeg_L) == (ps.soil.roof.ms,)
            @test size(exwater.ExWaterGroundImp_H) == (ps.soil.ground.ms,)
            @test size(exwater.ExWaterGroundImp_L) == (ps.soil.ground.ms,)
            @test size(exwater.ExWaterGroundBare_H) == (ps.soil.ground.ms,)
            @test size(exwater.ExWaterGroundBare_L) == (ps.soil.ground.ms,)
            @test size(exwater.ExWaterGroundVeg_H) == (ps.soil.ground.ms,)
            @test size(exwater.ExWaterGroundVeg_L) == (ps.soil.ground.ms,)
            @test size(exwater.ExWaterGroundTot_H) == (ps.soil.ground.ms,)
            @test size(exwater.ExWaterGroundTot_L) == (ps.soil.ground.ms,)

            # Test all fields are accessible, have correct dimensions and initialized to zero
            for field in fieldnames(ExWater)
                @test isa(getproperty(exwater, field), Array{FT,1})
                @test getproperty(exwater, field)[1] == 0
            end
        end

        @testset "ExWater matrix (N=2)" begin
            hours = 24
            exwater = initialize_exwater(FT, 2, ps.soil, hours)

            # Test structure
            @test exwater isa ExWater{FT,2}

            @test size(exwater.ExWaterRoofVeg_H) == (hours, ps.soil.roof.ms)
            @test size(exwater.ExWaterRoofVeg_L) == (hours, ps.soil.roof.ms)
            @test size(exwater.ExWaterGroundImp_H) == (hours, ps.soil.ground.ms)
            @test size(exwater.ExWaterGroundImp_L) == (hours, ps.soil.ground.ms)
            @test size(exwater.ExWaterGroundBare_H) == (hours, ps.soil.ground.ms)
            @test size(exwater.ExWaterGroundBare_L) == (hours, ps.soil.ground.ms)
            @test size(exwater.ExWaterGroundVeg_H) == (hours, ps.soil.ground.ms)
            @test size(exwater.ExWaterGroundVeg_L) == (hours, ps.soil.ground.ms)
            @test size(exwater.ExWaterGroundTot_H) == (hours, ps.soil.ground.ms)
            @test size(exwater.ExWaterGroundTot_L) == (hours, ps.soil.ground.ms)

            # Test all fields are accessible
            for field in fieldnames(ExWater)
                @test isa(getproperty(exwater, field), Matrix{FT})
            end
        end
    end

    @testset "SoilPotW" begin
        @testset "SoilPotW scalar (N=0)" begin
            soilpotw = initialize_soilpotw(FT, 0)

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

        @testset "SoilPotW vector (N=1)" begin
            hours = 24
            soilpotw = initialize_soilpotw(FT, 1, hours)

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
        @testset "CiCO2Leaf scalar (N=0)" begin
            cico2leaf = initialize_cico2leaf(FT, 0)

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

        @testset "CiCO2Leaf vector (N=1)" begin
            hours = 24
            initial_value = 300.0
            cico2leaf = initialize_cico2leaf(FT, 1, initial_value, hours)

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
    @testset "WaterFluxVariables 'scalar' (N=0)" begin
        water_flux_vars = initialize_water_flux_variables(FT, 0, ps.soil)

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

    @testset "WaterFluxVariables 'vector' (N=1)" begin
        # Mock the water_flux_variables initialization
        hours = 24
        initial_value = 300.0
        water_flux_vars = initialize_water_flux_variables(
            FT, 1, ps.soil, initial_value, hours
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
