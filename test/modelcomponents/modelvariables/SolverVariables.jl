using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    SolverVariables, initialize_solver_variables

FT = Float64

@testset "N = 0" begin
    sv = initialize_solver_variables(FT, 0)
    @test sv isa SolverVariables{FT,0,1}

    @test typeof(getfield(sv, :Success)) == Array{Bool,0}
    @test sv.Success == 0
end

@testset "N = 1" begin
    sv = initialize_solver_variables(FT, 1)
    @test sv isa SolverVariables{FT,1,2}
end
