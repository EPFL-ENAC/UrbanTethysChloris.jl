using Test
using UrbanTethysChloris.Water: water_impervious

FT = Float64

@testset "MATLAB" begin
    dth = 1.0
    E_imp = -2.281563103135239e-06
    In_imp_tm1 = 0.007787189104362
    In_max_imp = 0.25
    K_imp = 0.0
    Rain = 0.0
    row = 1000.0
    Runon_tm1 = 0.0

    q_runon_imp, In_imp, dIn_imp_dt, Lk_imp, WBalance_In_imp = water_impervious(
        Rain, Runon_tm1, E_imp, In_imp_tm1, dth, row, In_max_imp, K_imp
    )

    @test q_runon_imp ≈ 0
    @test In_imp ≈ 0.016000816275649
    @test dIn_imp_dt ≈ 0.008213627171287
    @test Lk_imp ≈ 0.0
    @test WBalance_In_imp ≈ 0.0 atol=eps(FT)
end
