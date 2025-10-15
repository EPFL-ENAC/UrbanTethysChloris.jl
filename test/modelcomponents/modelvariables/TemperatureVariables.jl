using Test
using TethysChlorisCore
using UrbanTethysChloris.ModelComponents.ModelVariables:
    TemperatureVariables, initialize_temperature_variables

FT = Float64

@testset "TemperatureVariables scalar (N=0)" begin
    temp_vars = initialize_temperature_variables(FT, 0)

    # Test structure
    @test temp_vars isa TemperatureVariables{FT,0}

    # Test field access for scalar case
    @test temp_vars.TRoofImp === 0.0
    @test temp_vars.UTCI === 0.0

    # Test all fields are accessible
    for field in fieldnames(TemperatureVariables)
        @test isa(getproperty(temp_vars, field), FT)
    end
end

@testset "TemperatureVariables vector (N=1)" begin
    hours = 24

    temp_vars = initialize_temperature_variables(FT, 1, hours)

    # Test structure
    @test temp_vars isa TemperatureVariables{FT,1}

    # Test field dimensions
    @test size(temp_vars.TRoofImp) == (hours,)
    @test size(temp_vars.UTCI) == (hours,)

    # Test all fields are accessible and have correct dimensions
    for field in fieldnames(TemperatureVariables)
        @test isa(getproperty(temp_vars, field), Array{FT,1})
        @test size(getproperty(temp_vars, field)) == (hours,)
    end
end
