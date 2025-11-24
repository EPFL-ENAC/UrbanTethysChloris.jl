using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Wind, LAITimeSeries, Resistance, EnvironmentalConditions

FT = Float64

@testset "Subsets" begin
    @testset "Wind" begin
        wind = Wind(FT)
        @test wind isa Wind{FT}
        for field in fieldnames(Wind)
            @test isa(getproperty(wind, field), FT)
            @test getproperty(wind, field) == 0
        end
    end

    @testset "LAITimeSeries" begin
        lai = LAITimeSeries(FT)
        @test lai isa LAITimeSeries{FT}
        for field in fieldnames(LAITimeSeries)
            @test isa(getproperty(lai, field), FT)
            @test getproperty(lai, field) == 0
        end
    end

    @testset "Resistance" begin
        res = Resistance(FT)
        @test res isa Resistance{FT}
        for field in fieldnames(Resistance)
            @test isa(getproperty(res, field), FT)
            @test getproperty(res, field) == 0
        end
    end
end

@testset "EnvironmentalConditions" begin
    env = EnvironmentalConditions(FT)
    @test env isa EnvironmentalConditions{FT}
    @test env.wind isa Wind{FT}
    @test env.LAI_time_series isa LAITimeSeries{FT}
    @test env.resistance isa Resistance{FT}
end
