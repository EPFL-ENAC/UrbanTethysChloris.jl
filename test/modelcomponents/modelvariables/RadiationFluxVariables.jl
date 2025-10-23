using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    AbsorbedRadiationFluxVariablesSubset,
    initialize_absorbed_radiation_flux_variables,
    DefaultRadiationFluxVariablesSubset,
    initialize_default_radiation_flux_variables,
    AlbedoOutput,
    initialize_albedo_output,
    RadiationFluxVariables,
    initialize_radiation_flux_variables,
    TimeSlice,
    TimeSeries

FT = Float64

@testset "Subsets" begin
    @testset "AbsorbedRadiationFluxVariablesSubset scalar (N=0)" begin
        swrabs = initialize_absorbed_radiation_flux_variables(FT, TimeSlice())

        # Test structure
        @test swrabs isa AbsorbedRadiationFluxVariablesSubset{FT,0}

        # Test all fields are accessible
        for field in fieldnames(AbsorbedRadiationFluxVariablesSubset)
            @test isa(getproperty(swrabs, field), FT)
            @test getproperty(swrabs, field) === 0.0
        end
    end

    @testset "AbsorbedRadiationFluxVariablesSubset vector (N=1)" begin
        hours = 24
        swrabs = initialize_absorbed_radiation_flux_variables(FT, TimeSeries(), hours)

        @test swrabs isa AbsorbedRadiationFluxVariablesSubset{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(AbsorbedRadiationFluxVariablesSubset)
            @test isa(getproperty(swrabs, field), Array{FT,1})
            @test size(getproperty(swrabs, field)) == (hours,)
            @test getproperty(swrabs, field)[1] == 0
        end
    end

    @testset "DefaultRadiationFluxVariablesSubset scalar (N=0)" begin
        swrin = initialize_default_radiation_flux_variables(FT, TimeSlice())

        # Test structure
        @test swrin isa DefaultRadiationFluxVariablesSubset{FT,0}

        # Test all fields are accessible
        for field in fieldnames(DefaultRadiationFluxVariablesSubset)
            @test isa(getproperty(swrin, field), FT)
            @test getproperty(swrin, field) === 0.0
        end
    end

    @testset "DefaultRadiationFluxVariablesSubset vector (N=1)" begin
        hours = 24
        swrin = initialize_default_radiation_flux_variables(FT, TimeSeries(), hours)

        @test swrin isa DefaultRadiationFluxVariablesSubset{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(DefaultRadiationFluxVariablesSubset)
            @test isa(getproperty(swrin, field), Array{FT,1})
            @test size(getproperty(swrin, field)) == (hours,)
            @test getproperty(swrin, field)[1] == 0
        end
    end

    @testset "AlbedoOutput scalar (N=0)" begin
        albedo_output = initialize_albedo_output(FT, TimeSlice())

        # Test structure
        @test albedo_output isa AlbedoOutput{FT,0}

        # Test all fields are accessible
        for field in fieldnames(AlbedoOutput)
            @test isa(getproperty(albedo_output, field), FT)
            @test getproperty(albedo_output, field) === 0.0
        end
    end

    @testset "AlbedoOutput vector (N=1)" begin
        hours = 24
        albedo_output = initialize_albedo_output(FT, TimeSeries(), hours)

        @test albedo_output isa AlbedoOutput{FT,1}

        # Test all fields are accessible, have correct dimensions and initialized to zero
        for field in fieldnames(AlbedoOutput)
            @test isa(getproperty(albedo_output, field), Array{FT,1})
            @test size(getproperty(albedo_output, field)) == (hours,)
            @test getproperty(albedo_output, field)[1] == 0
        end
    end
end

@testset "RadiationFluxVariables scalar (N=0)" begin
    radiation_vars = initialize_radiation_flux_variables(FT, TimeSlice())

    # Test structure
    @test radiation_vars isa RadiationFluxVariables{FT,0}

    # Test fields are properly initialized
    @test isa(radiation_vars.SWRabs, AbsorbedRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.SWRin, DefaultRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.SWRout, DefaultRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.SWREB, DefaultRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.LWRabs, DefaultRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.LWRin, DefaultRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.LWRout, DefaultRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.LWREB, DefaultRadiationFluxVariablesSubset{FT,0})
    @test isa(radiation_vars.AlbedoOutput, AlbedoOutput{FT,0})

    # Test some random fields
    @test radiation_vars.SWRabs.GroundImp === 0.0
    @test radiation_vars.LWRin.Tree === 0.0
    @test radiation_vars.AlbedoOutput.TotalUrban === 0.0
end

@testset "RadiationFluxVariables vector (N=1)" begin
    hours = 24
    radiation_vars = initialize_radiation_flux_variables(FT, TimeSeries(), hours)

    # Test structure
    @test radiation_vars isa RadiationFluxVariables{FT,1}

    # Test fields are properly initialized
    @test isa(radiation_vars.SWRabs, AbsorbedRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.SWRin, DefaultRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.SWRout, DefaultRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.SWREB, DefaultRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.LWRabs, DefaultRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.LWRin, DefaultRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.LWRout, DefaultRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.LWREB, DefaultRadiationFluxVariablesSubset{FT,1})
    @test isa(radiation_vars.AlbedoOutput, AlbedoOutput{FT,1})

    # Test dimensions of a few fields
    @test size(radiation_vars.SWRabs.GroundImp) == (hours,)
    @test size(radiation_vars.LWRin.Tree) == (hours,)
    @test size(radiation_vars.AlbedoOutput.TotalUrban) == (hours,)

    # Test initialization values
    @test radiation_vars.SWRabs.GroundImp[1] == 0.0
    @test radiation_vars.LWRin.Tree[1] == 0.0
    @test radiation_vars.AlbedoOutput.TotalUrban[1] == 0.0
end
