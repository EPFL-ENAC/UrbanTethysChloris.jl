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
    @testset "TimeSlice" begin
        wind = initialize_wind(FT, TimeSlice())
        @test wind isa Wind{FT,0}
        for field in fieldnames(Wind)
            @test isa(getproperty(wind, field), FT)
        end
    end

    @testset "TimeSeries" begin
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
    @testset "TimeSlice" begin
        lai = initialize_lai_time_series(FT, TimeSlice())
        @test lai isa LAITimeSeries{FT,0}
        for field in fieldnames(LAITimeSeries)
            @test isa(getproperty(lai, field), FT)
        end
    end

    @testset "TimeSeries" begin
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
    @testset "TimeSlice" begin
        res = initialize_resistance(FT, TimeSlice())
        @test res isa Resistance{FT,0}
        for field in fieldnames(Resistance)
            @test isa(getproperty(res, field), FT)
        end
    end

    @testset "TimeSeries" begin
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
    @testset "TimeSlice" begin
        env = initialize_environmental_conditions(FT, TimeSlice())
        @test env isa EnvironmentalConditions{FT,0}
        @test env.wind isa Wind{FT,0}
        @test env.LAI_time_series isa LAITimeSeries{FT,0}
        @test env.resistance isa Resistance{FT,0}
    end

    @testset "TimeSeries" begin
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

@testset "get/setindex" begin
    hours = 24
    env = initialize_environmental_conditions(FT, TimeSeries(), hours)

    LAI_R = FT(2.0)
    u_Hcan = FT(3.0)
    raRooftoAtm = FT(4.0)

    x = env[1]
    x.LAI_time_series.LAI_R = LAI_R
    x.wind.u_Hcan = u_Hcan
    x.resistance.raRooftoAtm = raRooftoAtm

    env[2] = x
    @test env.LAI_time_series.LAI_R[2] == LAI_R
    @test env.wind.u_Hcan[2] == u_Hcan
    @test env.resistance.raRooftoAtm[2] == raRooftoAtm
end
