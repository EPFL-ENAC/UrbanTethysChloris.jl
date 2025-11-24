using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    AbsorbedRadiationFluxVariablesSubset,
    DefaultRadiationFluxVariablesSubset,
    AlbedoOutput,
    RadiationFluxVariables

FT = Float64

@testset "Subsets" begin
    @testset "AbsorbedRadiationFluxVariablesSubset" begin
        swrabs = AbsorbedRadiationFluxVariablesSubset(FT)
        @test swrabs isa AbsorbedRadiationFluxVariablesSubset{FT}
        for field in fieldnames(AbsorbedRadiationFluxVariablesSubset)
            @test isa(getproperty(swrabs, field), FT)
            @test getproperty(swrabs, field) === 0.0
        end
    end

    @testset "DefaultRadiationFluxVariablesSubset" begin
        swrin = DefaultRadiationFluxVariablesSubset(FT)
        @test swrin isa DefaultRadiationFluxVariablesSubset{FT}
        for field in fieldnames(DefaultRadiationFluxVariablesSubset)
            @test isa(getproperty(swrin, field), FT)
            @test getproperty(swrin, field) === 0.0
        end
    end

    @testset "AlbedoOutput" begin
        albedo_output = AlbedoOutput(FT)
        @test albedo_output isa AlbedoOutput{FT}
        for field in fieldnames(AlbedoOutput)
            @test isa(getproperty(albedo_output, field), FT)
            @test getproperty(albedo_output, field) === 0.0
        end
    end
end

@testset "RadiationFluxVariables" begin
    radiation_vars = RadiationFluxVariables(FT)
    @test radiation_vars isa RadiationFluxVariables{FT}
    @test isa(radiation_vars.SWRabs, AbsorbedRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.SWRin, DefaultRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.SWRout, DefaultRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.SWREB, DefaultRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.LWRabs, DefaultRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.LWRin, DefaultRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.LWRout, DefaultRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.LWREB, DefaultRadiationFluxVariablesSubset{FT})
    @test isa(radiation_vars.AlbedoOutput, AlbedoOutput{FT})
    @test radiation_vars.SWRabs.GroundImp === 0.0
    @test radiation_vars.LWRin.Tree === 0.0
    @test radiation_vars.AlbedoOutput.TotalUrban === 0.0
end
