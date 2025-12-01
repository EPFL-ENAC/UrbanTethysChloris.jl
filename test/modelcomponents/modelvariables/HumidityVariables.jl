using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    Humidity, Results2m, HumidityVariables

FT = Float64

@testset "Subsets" begin
    @testset "Humidity" begin
        humidity = Humidity(FT)
        @test humidity isa Humidity{FT}
        @test humidity.CanyonSpecific === 0.0
        @test humidity.AtmSpecific === 0.0
        @test humidity.CanyonRelative === 0.0
        for field in fieldnames(Humidity)
            @test isa(getproperty(humidity, field), FT)
            @test getproperty(humidity, field) == 0
        end
    end

    @testset "Results2m" begin
        results = Results2m(FT)
        @test results isa Results2m{FT}
        @test results.e_T2m === 0.0
        for field in fieldnames(Results2m)
            @test isa(getproperty(results, field), FT)
            @test getproperty(results, field) == 0
        end
    end
end

@testset "HumidityVariables" begin
    humidity_vars = HumidityVariables(FT)
    @test humidity_vars isa HumidityVariables{FT}
    @test humidity_vars.Humidity isa Humidity{FT}
    @test humidity_vars.Results2m isa Results2m{FT}
end
