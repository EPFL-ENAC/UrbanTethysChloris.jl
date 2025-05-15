using Test
using Dates
using UrbanTethysChloris: set_sun_variables

Datam = DateTime(2010, 1, 1, 1, 0, 0)
DeltaGMT = 1.0
Lat = 47.380000000000003
Lon = 8.56
t_aft = 0.5
t_bef = 0.5

@testset "MATLAB" begin
    @testset "Basic case" begin
        h_S, zeta_S, T_sunrise, T_sunset = set_sun_variables(
            Datam, DeltaGMT, Lon, Lat, t_bef, t_aft
        )
        @test h_S ≈ -1.124694136568079
        @test zeta_S ≈ 0.321038765407712
        @test T_sunrise ≈ 7.830586347989140
        @test T_sunset ≈ 16.169413652010864
    end

    @testset "Night and negative longitude" begin
        h_S, zeta_S, T_sunrise, T_sunset = set_sun_variables(
            DateTime(2010, 1, 1, 23, 0, 0), DeltaGMT, -10.0, Lat, t_bef, t_aft
        )
        @test h_S ≈ -1.118980500706043
        @test zeta_S ≈ 5.910698991846242
        @test T_sunrise ≈ 7.830586347989140
        @test T_sunset ≈ 16.169413652010864
    end

    @testset "Morning" begin
        h_S, zeta_S, T_sunrise, T_sunset = set_sun_variables(
            DateTime(2010, 1, 1, 9, 0, 0), DeltaGMT, Lon, Lat, t_bef, t_aft
        )
        @test h_S ≈ 0.100212291676105
        @test zeta_S ≈ 2.333463084032588
        @test T_sunrise ≈ 7.830586347989140
        @test T_sunset ≈ 16.169413652010864
    end
end
