using Test
using StaticArrays
using UrbanTethysChloris.ModelComponents.ModelVariables:
    WBRoof, WBCanyonIndv, WBCanyonTot, EB, SolverVariables, EnergyBalanceVariables

FT = Float64

@testset "Subsets" begin
    @testset "WBRoof" begin
        wbroof = WBRoof(FT)
        @test wbroof isa WBRoof{FT}
        for field in fieldnames(WBRoof)
            @test isa(getproperty(wbroof, field), FT)
            @test getproperty(wbroof, field) == 0
        end
    end

    @testset "WBCanyonIndv" begin
        wbcanyonindv = WBCanyonIndv(FT)
        @test wbcanyonindv isa WBCanyonIndv{FT}
        for field in fieldnames(WBCanyonIndv)
            @test isa(getproperty(wbcanyonindv, field), FT)
            @test getproperty(wbcanyonindv, field) == 0
        end
    end

    @testset "WBCanyonTot" begin
        wbcanyontot = WBCanyonTot(FT)
        @test wbcanyontot isa WBCanyonTot{FT}
        for field in fieldnames(WBCanyonTot)
            @test isa(getproperty(wbcanyontot, field), FT)
            @test getproperty(wbcanyontot, field) == 0
        end
    end

    @testset "EB" begin
        eb = EB(FT)
        @test eb isa EB{FT}
        for field in fieldnames(EB)
            @test isa(getproperty(eb, field), FT)
            @test getproperty(eb, field) == 0
        end
    end

    @testset "SolverVariables" begin
        solver_vars = SolverVariables(FT)
        @test solver_vars isa SolverVariables{FT}
        @test solver_vars.Success === false
        @test solver_vars.ValuesEB isa Vector{FT}
        @test all(solver_vars.ValuesEB .== 0)
        @test solver_vars.Tsolver isa Vector{FT}
        @test all(solver_vars.Tsolver .== 0)
        @test solver_vars.YfunctionOutput isa Vector{FT}
        @test all(solver_vars.YfunctionOutput .== 0)
    end
end

@testset "EnergyBalanceVariables" begin
    eb_vars = EnergyBalanceVariables(FT)
    @test eb_vars isa EnergyBalanceVariables{FT}
    @test eb_vars.WBRoof isa WBRoof{FT}
    @test eb_vars.WBCanyonIndv isa WBCanyonIndv{FT}
    @test eb_vars.WBCanyonTot isa WBCanyonTot{FT}
    @test eb_vars.EB isa EB{FT}
    @test eb_vars.Solver isa SolverVariables{FT}
end
