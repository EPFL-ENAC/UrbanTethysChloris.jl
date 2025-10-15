using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    EnvironmentalConditions, initialize_environmental_conditions

FT = Float64

@testset "N = 0" begin
    ec = initialize_environmental_conditions(FT, 0)
    @test ec isa EnvironmentalConditions{FT,0}

    @test typeof(getfield(ec, :u_Hcan)) == Array{Float64,0}
    @test ec.u_Hcan == 0
end

@testset "N = 1" begin
    ec = initialize_environmental_conditions(FT, 1)
    @test ec isa EnvironmentalConditions{FT,1}
end
