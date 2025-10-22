using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Wind,
    initialize_wind,
    LAITimeSeries,
    initialize_lai_time_series,
    Resistance,
    initialize_resistance,
    EnvironmentalConditions,
    initialize_environmental_conditions

FT = Float64
hours = 24

@testset "Wind" begin
    @testset "N = 0" begin
        w = initialize_wind(FT, 0)
        @test w isa Wind{FT,0}

        for field in fieldnames(Wind)
            @test isa(getproperty(w, field), FT)
            @test getproperty(w, field) == 0
        end
    end

    @testset "N = 1" begin
        w = initialize_wind(FT, 1, hours)
        @test w isa Wind{FT,1}

        for field in fieldnames(Wind)
            @test isa(getproperty(w, field), Array{FT,1})
            @test size(getproperty(w, field)) == (hours,)
            @test getproperty(w, field)[1] == 0
        end
    end
end

@testset "LAITimeSeries" begin
    @testset "N = 0" begin
        LAI_ts = initialize_lai_time_series(FT, 0)
        @test LAI_ts isa LAITimeSeries{FT,0}

        for field in fieldnames(LAITimeSeries)
            @test isa(getproperty(LAI_ts, field), FT)
            @test getproperty(LAI_ts, field) == 0
        end
    end

    @testset "N = 1" begin
        LAI_ts = initialize_lai_time_series(FT, 1, hours)
        @test LAI_ts isa LAITimeSeries{FT,1}

        for field in fieldnames(LAITimeSeries)
            @test isa(getproperty(LAI_ts, field), Array{FT,1})
            @test size(getproperty(LAI_ts, field)) == (hours,)
            @test getproperty(LAI_ts, field)[1] == 0
        end
    end
end

@testset "Resistance" begin
    @testset "N = 0" begin
        RES = initialize_resistance(FT, 0)
        @test RES isa Resistance{FT,0}

        for field in fieldnames(Resistance)
            @test isa(getproperty(RES, field), FT)
            @test getproperty(RES, field) == 0
        end
    end

    @testset "N = 1" begin
        RES = initialize_resistance(FT, 1, hours)
        @test RES isa Resistance{FT,1}

        for field in fieldnames(Resistance)
            @test isa(getproperty(RES, field), Array{FT,1})
            @test size(getproperty(RES, field)) == (hours,)
            @test getproperty(RES, field)[1] == 0
        end
    end
end

@testset "EnvironmentalConditions" begin
    @testset "N = 0" begin
        ec = initialize_environmental_conditions(FT, 0)
        @test ec isa EnvironmentalConditions{FT,0}

        @test ec.wind isa Wind{FT,0}
        @test ec.LAI_time_series isa LAITimeSeries{FT,0}
        @test ec.resistance isa Resistance{FT,0}
    end

    @testset "N = 1" begin
        ec = initialize_environmental_conditions(FT, 1)
        @test ec isa EnvironmentalConditions{FT,1}

        @test ec.wind isa Wind{FT,1}
        @test ec.LAI_time_series isa LAITimeSeries{FT,1}
        @test ec.resistance isa Resistance{FT,1}
    end
end
