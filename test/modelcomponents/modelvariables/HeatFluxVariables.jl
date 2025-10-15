using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Hflux,
    initialize_hflux,
    LEflux,
    initialize_leflux,
    Gflux,
    initialize_gflux,
    dStorage,
    initialize_dstorage,
    Results2mEnergyFlux,
    initialize_results2m_energy_flux,
    HeatFluxVariables,
    initialize_heat_flux_variables

FT = Float64

@testset "Subsets" begin
    @testset "Hflux scalar (N=0)" begin
        hflux = initialize_hflux(FT, 0)

        # Test structure
        @test hflux isa Hflux{FT,0}

        # Test field access for scalar case
        @test hflux.HfluxRoofImp === 0.0
        @test hflux.HfluxCanyon === 0.0
        @test hflux.dS_H_air === 0.0

        # Test all fields are accessible
        for field in fieldnames(Hflux)
            @test isa(getproperty(hflux, field), FT)
        end
    end

    @testset "Hflux vector (N=1)" begin
        hours = 24
        hflux = initialize_hflux(FT, 1, hours)

        @test hflux isa Hflux{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Hflux)
            @test isa(getproperty(hflux, field), Array{FT,1})
            @test size(getproperty(hflux, field)) == (hours,)
            @test getproperty(hflux, field)[1] == 0
        end
    end

    @testset "LEflux scalar (N=0)" begin
        leflux = initialize_leflux(FT, 0)

        # Test structure
        @test leflux isa LEflux{FT,0}

        # Test field access for scalar case
        @test leflux.LEfluxRoofImp === 0.0
        @test leflux.LEfluxCanyon === 0.0
        @test leflux.dS_LE_air === 0.0

        # Test all fields are accessible
        for field in fieldnames(LEflux)
            @test isa(getproperty(leflux, field), FT)
        end
    end

    @testset "LEflux vector (N=1)" begin
        hours = 24
        leflux = initialize_leflux(FT, 1, hours)

        @test leflux isa LEflux{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(LEflux)
            @test isa(getproperty(leflux, field), Array{FT,1})
            @test size(getproperty(leflux, field)) == (hours,)
            @test getproperty(leflux, field)[1] == 0
        end
    end

    @testset "Gflux scalar (N=0)" begin
        gflux = initialize_gflux(FT, 0)

        # Test structure
        @test gflux isa Gflux{FT,0}

        # Test field access for scalar case
        @test gflux.G1RoofImp === 0.0
        @test gflux.G1Canyon === 0.0
        @test gflux.G2Urban === 0.0

        # Test all fields are accessible
        for field in fieldnames(Gflux)
            @test isa(getproperty(gflux, field), FT)
        end
    end

    @testset "Gflux vector (N=1)" begin
        hours = 24
        gflux = initialize_gflux(FT, 1, hours)

        @test gflux isa Gflux{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Gflux)
            @test isa(getproperty(gflux, field), Array{FT,1})
            @test size(getproperty(gflux, field)) == (hours,)
            @test getproperty(gflux, field)[1] == 0
        end
    end

    @testset "dStorage scalar (N=0)" begin
        dstorage = initialize_dstorage(FT, 0)

        # Test structure
        @test dstorage isa dStorage{FT,0}

        # Test field access for scalar case
        @test dstorage.dsRoofImp === 0.0
        @test dstorage.dsWallSun === 0.0
        @test dstorage.dsCanyonAir === 0.0

        # Test all fields are accessible
        for field in fieldnames(dStorage)
            @test isa(getproperty(dstorage, field), FT)
        end
    end

    @testset "dStorage vector (N=1)" begin
        hours = 24
        dstorage = initialize_dstorage(FT, 1, hours)

        @test dstorage isa dStorage{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(dStorage)
            @test isa(getproperty(dstorage, field), Array{FT,1})
            @test size(getproperty(dstorage, field)) == (hours,)
            @test getproperty(dstorage, field)[1] == 0
        end
    end

    @testset "Results2mEnergyFlux scalar (N=0)" begin
        results2m_energy_flux = initialize_results2m_energy_flux(FT, 0)

        # Test structure
        @test results2m_energy_flux isa Results2mEnergyFlux{FT,0}

        # Test field access for scalar case
        @test results2m_energy_flux.DHi === 0.0
        @test results2m_energy_flux.Hcan_2m === 0.0
        @test results2m_energy_flux.Ecan_2m === 0.0

        # Test all fields are accessible
        for field in fieldnames(Results2mEnergyFlux)
            @test isa(getproperty(results2m_energy_flux, field), FT)
        end
    end

    @testset "Results2mEnergyFlux vector (N=1)" begin
        hours = 24
        results2m_energy_flux = initialize_results2m_energy_flux(FT, 1, hours)

        @test results2m_energy_flux isa Results2mEnergyFlux{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(Results2mEnergyFlux)
            @test isa(getproperty(results2m_energy_flux, field), Array{FT,1})
            @test size(getproperty(results2m_energy_flux, field)) == (hours,)
            @test getproperty(results2m_energy_flux, field)[1] == 0
        end
    end
end

@testset "HeatFluxVariables scalar (N=0)" begin
    heat_flux_vars = initialize_heat_flux_variables(FT, 0)

    # Test structure
    @test heat_flux_vars isa HeatFluxVariables{FT,0}

    # Test fields are properly initialized
    @test isa(heat_flux_vars.Hflux, Hflux{FT,0})
    @test isa(heat_flux_vars.LEflux, LEflux{FT,0})
    @test isa(heat_flux_vars.Gflux, Gflux{FT,0})
    @test isa(heat_flux_vars.dStorage, dStorage{FT,0})
    @test isa(heat_flux_vars.Results2mEnergyFlux, Results2mEnergyFlux{FT,0})

    # Test some random fields
    @test heat_flux_vars.Hflux.HfluxRoofImp === 0.0
    @test heat_flux_vars.LEflux.LEfluxRoof === 0.0
    @test heat_flux_vars.Gflux.G1Urban === 0.0
    @test heat_flux_vars.dStorage.dsTree === 0.0
    @test heat_flux_vars.Results2mEnergyFlux.Hcan_2m === 0.0
end

@testset "HeatFluxVariables vector (N=1)" begin
    hours = 24
    heat_flux_vars = initialize_heat_flux_variables(FT, 1, hours)

    # Test structure
    @test heat_flux_vars isa HeatFluxVariables{FT,1}

    # Test fields are properly initialized
    @test isa(heat_flux_vars.Hflux, Hflux{FT,1})
    @test isa(heat_flux_vars.LEflux, LEflux{FT,1})
    @test isa(heat_flux_vars.Gflux, Gflux{FT,1})
    @test isa(heat_flux_vars.dStorage, dStorage{FT,1})
    @test isa(heat_flux_vars.Results2mEnergyFlux, Results2mEnergyFlux{FT,1})

    # Test dimensions of a few fields
    @test size(heat_flux_vars.Hflux.HfluxRoofImp) == (hours,)
    @test size(heat_flux_vars.LEflux.LEfluxRoof) == (hours,)
    @test size(heat_flux_vars.Gflux.G1Urban) == (hours,)
    @test size(heat_flux_vars.dStorage.dsTree) == (hours,)
    @test size(heat_flux_vars.Results2mEnergyFlux.Hcan_2m) == (hours,)

    # Test initialization values
    @test heat_flux_vars.Hflux.HfluxRoofImp[1] == 0.0
    @test heat_flux_vars.LEflux.LEfluxRoof[1] == 0.0
    @test heat_flux_vars.Gflux.G1Urban[1] == 0.0
    @test heat_flux_vars.dStorage.dsTree[1] == 0.0
    @test heat_flux_vars.Results2mEnergyFlux.Hcan_2m[1] == 0.0
end
