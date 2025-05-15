using Test
using UrbanTethysChloris: incoming_longwave

Ta = 25.0
ea = 1000.0

@testset "MATLAB" begin
    @testset "Cloudiness" begin
        Latm = incoming_longwave(Ta, ea, 0.5)
        @test Latm ≈ 378.5560335493908
    end

    @testset "Incoming Longwave Radiation" begin
        Latm = incoming_longwave(Ta, ea, 1000.0)
        @test Latm ≈ 1000.0
    end
end
