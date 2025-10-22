using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Humidity,
    initialize_humidity,
    Results2m,
    initialize_results2m,
    HumidityVariables,
    initialize_humidity_variables

FT = Float64

@testset "Humidity" begin
    @testset "Humidity scalar (N=0)" begin
        atm_specific = FT(0.01)  # Initial specific humidity
        humidity = initialize_humidity(FT, 0, 1, atm_specific)

        # Test structure
        @test humidity isa Humidity{FT,0}

        # Test field access for scalar case
        @test humidity.CanyonSpecific === 0.0
        @test humidity.AtmSpecific === 0.0
        @test humidity.CanyonRelative === 0.0

        # Test all fields are accessible
        for field in fieldnames(Humidity)
            @test isa(getproperty(humidity, field), FT)
        end
    end

    @testset "Humidity vector (N=1)" begin
        hours = 24
        atm_specific = FT(0.01)  # Initial specific humidity

        humidity = initialize_humidity(FT, 1, hours, atm_specific)

        # Test structure
        @test humidity isa Humidity{FT,1}

        # Test field dimensions
        @test size(humidity.CanyonRelative) == (hours,)
        @test size(humidity.AtmSpecific) == (hours,)

        # Test initialization
        @test humidity.CanyonSpecific[1] == atm_specific
        @test humidity.AtmSpecific[1] == 0

        # Test all fields are accessible and have correct dimensions
        for field in fieldnames(Humidity)
            @test isa(getproperty(humidity, field), Array{FT,1})
            @test size(getproperty(humidity, field)) == (hours,)
        end
    end
end

@testset "Results2m" begin
    @testset "Results2m scalar (N=0)" begin
        hours = 24
        qatm = FT(0.01)   # Initial specific humidity

        results = initialize_results2m(FT, 0)

        # Test structure
        @test results isa Results2m{FT,0}

        # Test field access for scalar case
        @test results.e_T2m === 0.0

        # Test all fields are accessible
        for field in fieldnames(Results2m)
            @test isa(getproperty(results, field), FT)
        end
    end

    @testset "Results2m vector (N=1)" begin
        hours = 24

        results = initialize_results2m(FT, 1, hours)

        # Test structure
        @test results isa Results2m{FT,1}

        # Test field dimensions
        @test size(results.T2m) == (hours,)
        @test size(results.q2m) == (hours,)

        # Test initialization
        @test all(results.T2m .== 0)

        # Test all fields are accessible and have correct dimensions
        for field in fieldnames(Results2m)
            @test isa(getproperty(results, field), Array{FT,1})
            @test size(getproperty(results, field)) == (hours,)
        end
    end
end

@testset "HumidityVariables" begin
    @testset "HumidityVariables scalar (N=0)" begin
        humidity_vars = initialize_humidity_variables(FT, 0)

        # Test structure
        @test humidity_vars isa HumidityVariables{FT,0}
        @test humidity_vars.Humidity isa Humidity{FT,0}
        @test humidity_vars.Results2m isa Results2m{FT,0}
    end

    @testset "HumidityVariables vector (N=1)" begin
        hours = 24
        q_atm = 0.005

        humidity_vars = initialize_humidity_variables(FT, 1, hours, q_atm)

        # Test structure
        @test humidity_vars isa HumidityVariables{FT,1}
        @test humidity_vars.Humidity isa Humidity{FT,1}
        @test humidity_vars.Results2m isa Results2m{FT,1}

        # Test field dimensions
        @test size(humidity_vars.Humidity.CanyonRelative) == (hours,)
        @test size(humidity_vars.Results2m.T2m) == (hours,)

        # Test initial condition for CanyonSpecific
        @test humidity_vars.Humidity.CanyonSpecific[1] == q_atm
    end
end
