using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Hflux, LEflux, Gflux, dStorage, Results2mEnergyFlux, HeatFluxVariables

FT = Float64

@testset "Subsets" begin
    @testset "Hflux" begin
        hflux = Hflux(FT)
        @test hflux isa Hflux{FT}
        @test hflux.HfluxRoofImp === 0.0
        @test hflux.HfluxCanyon === 0.0
        @test hflux.dS_H_air === 0.0
        for field in fieldnames(Hflux)
            @test isa(getproperty(hflux, field), FT)
        end
    end

    @testset "LEflux" begin
        leflux = LEflux(FT)
        @test leflux isa LEflux{FT}
        @test leflux.LEfluxRoofImp === 0.0
        @test leflux.LEfluxCanyon === 0.0
        @test leflux.dS_LE_air === 0.0
        for field in fieldnames(LEflux)
            @test isa(getproperty(leflux, field), FT)
        end
    end

    @testset "Gflux" begin
        gflux = Gflux(FT)
        @test gflux isa Gflux{FT}
        @test gflux.G1RoofImp === 0.0
        @test gflux.G1Canyon === 0.0
        @test gflux.G2Urban === 0.0
        for field in fieldnames(Gflux)
            @test isa(getproperty(gflux, field), FT)
        end
    end

    @testset "dStorage" begin
        dstorage = dStorage(FT)
        @test dstorage isa dStorage{FT}
        @test dstorage.dsRoofImp === 0.0
        @test dstorage.dsWallSun === 0.0
        @test dstorage.dsCanyonAir === 0.0
        for field in fieldnames(dStorage)
            @test isa(getproperty(dstorage, field), FT)
        end
    end

    @testset "Results2mEnergyFlux" begin
        results2m_energy_flux = Results2mEnergyFlux(FT)
        @test results2m_energy_flux isa Results2mEnergyFlux{FT}
        @test results2m_energy_flux.DHi === 0.0
        @test results2m_energy_flux.Hcan_2m === 0.0
        @test results2m_energy_flux.Ecan_2m === 0.0
        for field in fieldnames(Results2mEnergyFlux)
            @test isa(getproperty(results2m_energy_flux, field), FT)
        end
    end
end

@testset "HeatFluxVariables" begin
    heat_flux_vars = HeatFluxVariables(FT)
    @test heat_flux_vars isa HeatFluxVariables{FT}
    @test isa(heat_flux_vars.Hflux, Hflux{FT})
    @test isa(heat_flux_vars.LEflux, LEflux{FT})
    @test isa(heat_flux_vars.Gflux, Gflux{FT})
    @test isa(heat_flux_vars.dStorage, dStorage{FT})
    @test isa(heat_flux_vars.Results2mEnergyFlux, Results2mEnergyFlux{FT})
    @test heat_flux_vars.Hflux.HfluxRoofImp === 0.0
    @test heat_flux_vars.LEflux.LEfluxRoof === 0.0
    @test heat_flux_vars.Gflux.G1Urban === 0.0
    @test heat_flux_vars.dStorage.dsTree === 0.0
    @test heat_flux_vars.Results2mEnergyFlux.Hcan_2m === 0.0
end
