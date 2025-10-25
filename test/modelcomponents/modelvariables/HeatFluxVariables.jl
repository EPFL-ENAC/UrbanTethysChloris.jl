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
    initialize_heat_flux_variables,
    TimeSlice,
    TimeSeries

FT = Float64

@testset "Subsets" begin
    @testset "Hflux" begin
        @testset "TimeSlice" begin
            hflux = initialize_hflux(FT, TimeSlice())
            @test hflux isa Hflux{FT,0}
            @test hflux.HfluxRoofImp === 0.0
            @test hflux.HfluxCanyon === 0.0
            @test hflux.dS_H_air === 0.0
            for field in fieldnames(Hflux)
                @test isa(getproperty(hflux, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            hflux = initialize_hflux(FT, TimeSeries(), hours)
            @test hflux isa Hflux{FT,1}
            for field in fieldnames(Hflux)
                @test isa(getproperty(hflux, field), Array{FT,1})
                @test size(getproperty(hflux, field)) == (hours,)
                @test getproperty(hflux, field)[1] == 0
            end
        end
    end

    @testset "LEflux" begin
        @testset "TimeSlice" begin
            leflux = initialize_leflux(FT, TimeSlice())
            @test leflux isa LEflux{FT,0}
            @test leflux.LEfluxRoofImp === 0.0
            @test leflux.LEfluxCanyon === 0.0
            @test leflux.dS_LE_air === 0.0
            for field in fieldnames(LEflux)
                @test isa(getproperty(leflux, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            leflux = initialize_leflux(FT, TimeSeries(), hours)
            @test leflux isa LEflux{FT,1}
            for field in fieldnames(LEflux)
                @test isa(getproperty(leflux, field), Array{FT,1})
                @test size(getproperty(leflux, field)) == (hours,)
                @test getproperty(leflux, field)[1] == 0
            end
        end
    end

    @testset "Gflux" begin
        @testset "TimeSlice" begin
            gflux = initialize_gflux(FT, TimeSlice())
            @test gflux isa Gflux{FT,0}
            @test gflux.G1RoofImp === 0.0
            @test gflux.G1Canyon === 0.0
            @test gflux.G2Urban === 0.0
            for field in fieldnames(Gflux)
                @test isa(getproperty(gflux, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            gflux = initialize_gflux(FT, TimeSeries(), hours)
            @test gflux isa Gflux{FT,1}
            for field in fieldnames(Gflux)
                @test isa(getproperty(gflux, field), Array{FT,1})
                @test size(getproperty(gflux, field)) == (hours,)
                @test getproperty(gflux, field)[1] == 0
            end
        end
    end

    @testset "dStorage" begin
        @testset "TimeSlice" begin
            dstorage = initialize_dstorage(FT, TimeSlice())
            @test dstorage isa dStorage{FT,0}
            @test dstorage.dsRoofImp === 0.0
            @test dstorage.dsWallSun === 0.0
            @test dstorage.dsCanyonAir === 0.0
            for field in fieldnames(dStorage)
                @test isa(getproperty(dstorage, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            dstorage = initialize_dstorage(FT, TimeSeries(), hours)
            @test dstorage isa dStorage{FT,1}
            for field in fieldnames(dStorage)
                @test isa(getproperty(dstorage, field), Array{FT,1})
                @test size(getproperty(dstorage, field)) == (hours,)
                @test getproperty(dstorage, field)[1] == 0
            end
        end
    end

    @testset "Results2mEnergyFlux" begin
        @testset "TimeSlice" begin
            results2m_energy_flux = initialize_results2m_energy_flux(FT, TimeSlice())
            @test results2m_energy_flux isa Results2mEnergyFlux{FT,0}
            @test results2m_energy_flux.DHi === 0.0
            @test results2m_energy_flux.Hcan_2m === 0.0
            @test results2m_energy_flux.Ecan_2m === 0.0
            for field in fieldnames(Results2mEnergyFlux)
                @test isa(getproperty(results2m_energy_flux, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            results2m_energy_flux = initialize_results2m_energy_flux(
                FT, TimeSeries(), hours
            )
            @test results2m_energy_flux isa Results2mEnergyFlux{FT,1}
            for field in fieldnames(Results2mEnergyFlux)
                @test isa(getproperty(results2m_energy_flux, field), Array{FT,1})
                @test size(getproperty(results2m_energy_flux, field)) == (hours,)
                @test getproperty(results2m_energy_flux, field)[1] == 0
            end
        end
    end
end

@testset "HeatFluxVariables" begin
    @testset "TimeSlice" begin
        heat_flux_vars = initialize_heat_flux_variables(FT, TimeSlice())
        @test heat_flux_vars isa HeatFluxVariables{FT,0}
        @test isa(heat_flux_vars.Hflux, Hflux{FT,0})
        @test isa(heat_flux_vars.LEflux, LEflux{FT,0})
        @test isa(heat_flux_vars.Gflux, Gflux{FT,0})
        @test isa(heat_flux_vars.dStorage, dStorage{FT,0})
        @test isa(heat_flux_vars.Results2mEnergyFlux, Results2mEnergyFlux{FT,0})
        @test heat_flux_vars.Hflux.HfluxRoofImp === 0.0
        @test heat_flux_vars.LEflux.LEfluxRoof === 0.0
        @test heat_flux_vars.Gflux.G1Urban === 0.0
        @test heat_flux_vars.dStorage.dsTree === 0.0
        @test heat_flux_vars.Results2mEnergyFlux.Hcan_2m === 0.0
    end

    @testset "TimeSeries" begin
        hours = 24
        heat_flux_vars = initialize_heat_flux_variables(FT, TimeSeries(), hours)
        @test heat_flux_vars isa HeatFluxVariables{FT,1}
        @test isa(heat_flux_vars.Hflux, Hflux{FT,1})
        @test isa(heat_flux_vars.LEflux, LEflux{FT,1})
        @test isa(heat_flux_vars.Gflux, Gflux{FT,1})
        @test isa(heat_flux_vars.dStorage, dStorage{FT,1})
        @test isa(heat_flux_vars.Results2mEnergyFlux, Results2mEnergyFlux{FT,1})
        @test size(heat_flux_vars.Hflux.HfluxRoofImp) == (hours,)
        @test size(heat_flux_vars.LEflux.LEfluxRoof) == (hours,)
        @test size(heat_flux_vars.Gflux.G1Urban) == (hours,)
        @test size(heat_flux_vars.dStorage.dsTree) == (hours,)
        @test size(heat_flux_vars.Results2mEnergyFlux.Hcan_2m) == (hours,)
        @test heat_flux_vars.Hflux.HfluxRoofImp[1] == 0.0
        @test heat_flux_vars.LEflux.LEfluxRoof[1] == 0.0
        @test heat_flux_vars.Gflux.G1Urban[1] == 0.0
        @test heat_flux_vars.dStorage.dsTree[1] == 0.0
        @test heat_flux_vars.Results2mEnergyFlux.Hcan_2m[1] == 0.0
    end
end
