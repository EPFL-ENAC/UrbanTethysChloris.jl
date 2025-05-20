using Test
using UrbanTethysChloris.Water: water_vegetation

FT = Float64

@testset "MATLAB" begin
    @testset "Zurich" begin
        dth = 1.0
        E_veg = 0.0
        In_veg_tm1 = 0.0
        LAI = 2.5
        Rain = 0.0
        row = 1000.0
        SAI = 0.001
        Sp_In = 0.2

        q_runon_veg, In_veg, dIn_veg_dt, WBalance_In_veg = water_vegetation(
            Rain, E_veg, In_veg_tm1, Sp_In, LAI, SAI, row, dth
        )

        @test q_runon_veg == 0
        @test In_veg == 0
        @test dIn_veg_dt == 0
        @test WBalance_In_veg == 0
    end

    @testset "Random" begin
        dth = 1.0
        E_veg = -2.281563103135239e-06
        In_veg_tm1 = 0.007787189104362
        LAI = 2.0
        Rain = 0.0
        row = 1000.0
        SAI = 0.5
        Sp_In = 0.2

        q_runon_veg, In_veg, dIn_veg_dt, WBalance_In_veg = water_vegetation(
            Rain, E_veg, In_veg_tm1, Sp_In, LAI, SAI, row, dth
        )

        @test q_runon_veg ≈ 0.010009629527896
        @test In_veg ≈ 0.005991186747753
        @test dIn_veg_dt ≈ -0.001796002356609
        @test WBalance_In_veg ≈ 0.0 atol=eps(FT)
    end
end
