using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Humidity,
    initialize_humidity,
    Results2m,
    initialize_results2m,
    HumidityVariables,
    initialize_humidity_variables,
    TimeSlice,
    TimeSeries

FT = Float64

@testset "Humidity" begin
    @testset "TimeSlice" begin
        humidity = initialize_humidity(FT, TimeSlice())
        @test humidity isa Humidity{FT,0}
        @test humidity.CanyonSpecific === 0.0
        @test humidity.AtmSpecific === 0.0
        @test humidity.CanyonRelative === 0.0
        for field in fieldnames(Humidity)
            @test isa(getproperty(humidity, field), FT)
        end
    end

    @testset "TimeSeries" begin
        hours = 24
        atm_specific = FT(0.01)
        humidity = initialize_humidity(FT, TimeSeries(), hours, atm_specific)
        @test humidity isa Humidity{FT,1}
        @test size(humidity.CanyonRelative) == (hours,)
        @test size(humidity.AtmSpecific) == (hours,)
        @test humidity.CanyonSpecific[1] == atm_specific
        @test humidity.AtmSpecific[1] == 0
        for field in fieldnames(Humidity)
            @test isa(getproperty(humidity, field), Array{FT,1})
            @test size(getproperty(humidity, field)) == (hours,)
        end
    end
end

@testset "Results2m" begin
    @testset "TimeSlice" begin
        results = initialize_results2m(FT, TimeSlice())
        @test results isa Results2m{FT,0}
        @test results.e_T2m === 0.0
        for field in fieldnames(Results2m)
            @test isa(getproperty(results, field), FT)
        end
    end

    @testset "TimeSeries" begin
        hours = 24
        results = initialize_results2m(FT, TimeSeries(), hours)
        @test results isa Results2m{FT,1}
        @test size(results.T2m) == (hours,)
        @test size(results.q2m) == (hours,)
        @test all(results.T2m .== 0)
        for field in fieldnames(Results2m)
            @test isa(getproperty(results, field), Array{FT,1})
            @test size(getproperty(results, field)) == (hours,)
        end
    end
end

@testset "HumidityVariables" begin
    @testset "TimeSlice" begin
        humidity_vars = initialize_humidity_variables(FT, TimeSlice())
        @test humidity_vars isa HumidityVariables{FT,0}
        @test humidity_vars.Humidity isa Humidity{FT,0}
        @test humidity_vars.Results2m isa Results2m{FT,0}
    end

    @testset "TimeSeries" begin
        hours = 24
        q_atm = 0.005
        humidity_vars = initialize_humidity_variables(FT, TimeSeries(), hours, q_atm)
        @test humidity_vars isa HumidityVariables{FT,1}
        @test humidity_vars.Humidity isa Humidity{FT,1}
        @test humidity_vars.Results2m isa Results2m{FT,1}
        @test size(humidity_vars.Humidity.CanyonRelative) == (hours,)
        @test size(humidity_vars.Results2m.T2m) == (hours,)
        @test humidity_vars.Humidity.CanyonSpecific[1] == q_atm
    end
end
