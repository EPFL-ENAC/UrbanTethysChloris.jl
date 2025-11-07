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
    @testset "AbsorbedRadiationFluxVariablesSubset" begin
        @testset "TimeSlice" begin
            swrabs = initialize_absorbed_radiation_flux_variables(FT, TimeSlice())
            @test swrabs isa AbsorbedRadiationFluxVariablesSubset{FT,0}
            for field in fieldnames(AbsorbedRadiationFluxVariablesSubset)
                @test isa(getproperty(swrabs, field), FT)
                @test getproperty(swrabs, field) === 0.0
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            swrabs = initialize_absorbed_radiation_flux_variables(FT, TimeSeries(), hours)
            @test swrabs isa AbsorbedRadiationFluxVariablesSubset{FT,1}
            for field in fieldnames(AbsorbedRadiationFluxVariablesSubset)
                @test isa(getproperty(swrabs, field), Array{FT,1})
                @test size(getproperty(swrabs, field)) == (hours,)
                @test getproperty(swrabs, field)[1] == 0
            end
        end
    end

    @testset "DefaultRadiationFluxVariablesSubset" begin
        @testset "TimeSlice" begin
            swrin = initialize_default_radiation_flux_variables(FT, TimeSlice())
            @test swrin isa DefaultRadiationFluxVariablesSubset{FT,0}
            for field in fieldnames(DefaultRadiationFluxVariablesSubset)
                @test isa(getproperty(swrin, field), FT)
                @test getproperty(swrin, field) === 0.0
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            swrin = initialize_default_radiation_flux_variables(FT, TimeSeries(), hours)
            @test swrin isa DefaultRadiationFluxVariablesSubset{FT,1}
            for field in fieldnames(DefaultRadiationFluxVariablesSubset)
                @test isa(getproperty(swrin, field), Array{FT,1})
                @test size(getproperty(swrin, field)) == (hours,)
                @test getproperty(swrin, field)[1] == 0
            end
        end
    end

    @testset "AlbedoOutput" begin
        @testset "TimeSlice" begin
            albedo_output = initialize_albedo_output(FT, TimeSlice())
            @test albedo_output isa AlbedoOutput{FT,0}
            for field in fieldnames(AlbedoOutput)
                @test isa(getproperty(albedo_output, field), FT)
                @test getproperty(albedo_output, field) === 0.0
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            albedo_output = initialize_albedo_output(FT, TimeSeries(), hours)
            @test albedo_output isa AlbedoOutput{FT,1}
            for field in fieldnames(AlbedoOutput)
                @test isa(getproperty(albedo_output, field), Array{FT,1})
                @test size(getproperty(albedo_output, field)) == (hours,)
                @test getproperty(albedo_output, field)[1] == 0
            end
        end
    end
end

@testset "RadiationFluxVariables" begin
    @testset "TimeSlice" begin
        radiation_vars = initialize_radiation_flux_variables(FT, TimeSlice())
        @test radiation_vars isa RadiationFluxVariables{FT,0}
        @test isa(radiation_vars.SWRabs, AbsorbedRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.SWRin, DefaultRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.SWRout, DefaultRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.SWREB, DefaultRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.LWRabs, DefaultRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.LWRin, DefaultRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.LWRout, DefaultRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.LWREB, DefaultRadiationFluxVariablesSubset{FT,0})
        @test isa(radiation_vars.AlbedoOutput, AlbedoOutput{FT,0})
        @test radiation_vars.SWRabs.GroundImp === 0.0
        @test radiation_vars.LWRin.Tree === 0.0
        @test radiation_vars.AlbedoOutput.TotalUrban === 0.0
    end

    @testset "TimeSeries" begin
        hours = 24
        radiation_vars = initialize_radiation_flux_variables(FT, TimeSeries(), hours)
        @test radiation_vars isa RadiationFluxVariables{FT,1}
        @test isa(radiation_vars.SWRabs, AbsorbedRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.SWRin, DefaultRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.SWRout, DefaultRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.SWREB, DefaultRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.LWRabs, DefaultRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.LWRin, DefaultRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.LWRout, DefaultRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.LWREB, DefaultRadiationFluxVariablesSubset{FT,1})
        @test isa(radiation_vars.AlbedoOutput, AlbedoOutput{FT,1})
        @test size(radiation_vars.SWRabs.GroundImp) == (hours,)
        @test size(radiation_vars.LWRin.Tree) == (hours,)
        @test size(radiation_vars.AlbedoOutput.TotalUrban) == (hours,)
        @test radiation_vars.SWRabs.GroundImp[1] == 0.0
        @test radiation_vars.LWRin.Tree[1] == 0.0
        @test radiation_vars.AlbedoOutput.TotalUrban[1] == 0.0
    end
end

@testset "get/setindex" begin
    hours = 24
    radiation_vars = initialize_radiation_flux_variables(FT, TimeSeries(), hours)

    SWRabs_GroundImp = FT(1.0)
    SWRin_Tree = FT(2.0)
    SWRout_GroundBare = FT(3.0)
    SWREB_TotalRoof = FT(4.0)
    LWRabs_GroundVeg = FT(5.0)
    LWRin_Tree = FT(6.0)
    LWRout_TotalGround = FT(7.0)
    LWREB_TotalCanyon = FT(8.0)
    AlbedoOutput_Roof = FT(0.5)

    x = radiation_vars[1]
    x.SWRabs.GroundImp = SWRabs_GroundImp
    x.SWRin.Tree = SWRin_Tree
    x.SWRout.GroundBare = SWRout_GroundBare
    x.SWREB.TotalRoof = SWREB_TotalRoof
    x.LWRabs.GroundVeg = LWRabs_GroundVeg
    x.LWRin.Tree = LWRin_Tree
    x.LWRout.TotalGround = LWRout_TotalGround
    x.LWREB.TotalCanyon = LWREB_TotalCanyon
    x.AlbedoOutput.Roof = AlbedoOutput_Roof

    radiation_vars[2] = x
    @test radiation_vars.SWRabs.GroundImp[2] == SWRabs_GroundImp
    @test radiation_vars.SWRin.Tree[2] == SWRin_Tree
    @test radiation_vars.SWRout.GroundBare[2] == SWRout_GroundBare
    @test radiation_vars.SWREB.TotalRoof[2] == SWREB_TotalRoof
    @test radiation_vars.LWRabs.GroundVeg[2] == LWRabs_GroundVeg
    @test radiation_vars.LWRin.Tree[2] == LWRin_Tree
    @test radiation_vars.LWRout.TotalGround[2] == LWRout_TotalGround
    @test radiation_vars.LWREB.TotalCanyon[2] == LWREB_TotalCanyon
    @test radiation_vars.AlbedoOutput.Roof[2] == AlbedoOutput_Roof
end
