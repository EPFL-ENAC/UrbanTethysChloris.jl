using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Wind,
    initialize_wind,
    LAITimeSeries,
    initialize_lai_time_series,
    Resistance,
    initialize_resistance,
    EnvironmentalConditions,
    initialize_environmental_conditions,
    TimeSlice,
    TimeSeries

FT = Float64

@testset "Wind" begin
    @testset "Wind scalar (N=0)" begin
        wind = initialize_wind(FT, TimeSlice())
        @test wind isa Wind{FT,0}
        for field in fieldnames(Wind)
            @test isa(getproperty(wind, field), FT)
        end
    end

    @testset "Wind vector (N=1)" begin
        hours = 24
        wind = initialize_wind(FT, TimeSeries(), hours)
        @test wind isa Wind{FT,1}
        for field in fieldnames(Wind)
            @test isa(getproperty(wind, field), Array{FT,1})
            @test size(getproperty(wind, field)) == (hours,)
        end
    end
end

@testset "LAITimeSeries" begin
    @testset "LAITimeSeries scalar (N=0)" begin
        lai = initialize_lai_time_series(FT, TimeSlice())
        @test lai isa LAITimeSeries{FT,0}
        for field in fieldnames(LAITimeSeries)
            @test isa(getproperty(lai, field), FT)
        end
    end

    @testset "LAITimeSeries vector (N=1)" begin
        hours = 24
        lai = initialize_lai_time_series(FT, TimeSeries(), hours)
        @test lai isa LAITimeSeries{FT,1}
        for field in fieldnames(LAITimeSeries)
            @test isa(getproperty(lai, field), Array{FT,1})
            @test size(getproperty(lai, field)) == (hours,)
        end
    end
end

@testset "Resistance" begin
    @testset "Resistance scalar (N=0)" begin
        res = initialize_resistance(FT, TimeSlice())
        @test res isa Resistance{FT,0}
        for field in fieldnames(Resistance)
            @test isa(getproperty(res, field), FT)
        end
    end

    @testset "Resistance vector (N=1)" begin
        hours = 24
        res = initialize_resistance(FT, TimeSeries(), hours)
        @test res isa Resistance{FT,1}
        for field in fieldnames(Resistance)
            @test isa(getproperty(res, field), Array{FT,1})
            @test size(getproperty(res, field)) == (hours,)
        end
    end
end

@testset "EnvironmentalConditions" begin
    @testset "EnvironmentalConditions scalar (N=0)" begin
        env = initialize_environmental_conditions(FT, TimeSlice())
        @test env isa EnvironmentalConditions{FT,0}
        @test env.wind isa Wind{FT,0}
        @test env.LAI_time_series isa LAITimeSeries{FT,0}
        @test env.resistance isa Resistance{FT,0}
    end

    @testset "EnvironmentalConditions vector (N=1)" begin
        hours = 24
        env = initialize_environmental_conditions(FT, TimeSeries(), hours)
        @test env isa EnvironmentalConditions{FT,1}
        @test env.wind isa Wind{FT,1}
        @test env.LAI_time_series isa LAITimeSeries{FT,1}
        @test env.resistance isa Resistance{FT,1}
        @test size(env.wind.u_Hcan) == (hours,)
        @test size(env.LAI_time_series.LAI_R) == (hours,)
        @test size(env.resistance.raRooftoAtm) == (hours,)
    end
end
