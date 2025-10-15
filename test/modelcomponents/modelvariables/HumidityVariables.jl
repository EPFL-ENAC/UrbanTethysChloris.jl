using Test
using TethysChlorisCore
using UrbanTethysChloris.ModelComponents.ModelVariables:
    HumidityVariables, initialize_humidity_variables

FT = Float64

@testset "HumidityVariables scalar (N=0)" begin
    humidity_vars = initialize_humidity_variables(FT, 0)

    # Test structure
    @test humidity_vars isa HumidityVariables{FT,0}

    # Test field access for scalar case
    @test humidity_vars.CanyonRelative === 0.0
    @test humidity_vars.AtmSpecific === 0.0
    @test humidity_vars.T2m === 0.0
    @test humidity_vars.RH_Tcan === 0.0

    # Test all fields are accessible
    for field in fieldnames(HumidityVariables)
        @test isa(getproperty(humidity_vars, field), FT)
    end
end

@testset "HumidityVariables vector (N=1)" begin
    hours = 24
    q_atm = 0.005

    humidity_vars = initialize_humidity_variables(FT, 1, hours, q_atm)

    # Test structure
    @test humidity_vars isa HumidityVariables{FT,1}

    # Test field dimensions
    @test size(humidity_vars.CanyonRelative) == (hours,)
    @test size(humidity_vars.T2m) == (hours,)

    # Test all fields are accessible and have correct dimensions
    for field in fieldnames(HumidityVariables)
        @test isa(getproperty(humidity_vars, field), Array{FT,1})
        @test size(getproperty(humidity_vars, field)) == (hours,)
    end

    # Test initial condition for CanyonSpecific
    @test humidity_vars.CanyonSpecific[1] == q_atm
end
