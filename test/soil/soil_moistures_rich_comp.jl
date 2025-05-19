using Test
using UrbanTethysChloris.Soil: soil_moistures_rich_comp

FT = Float64

@testset "MATLAB" begin
    alpVG = fill(-9.634329919190001e-05, 4)
    aR = 1.0
    aT = 1000.0
    cosalp = 1.0
    dz = [10.0, 10.0, 30.0, 50.0]
    Dz = [5.0, 10.0, 20.0, 40.0]
    EG = [2.068630444258776e-06, 0.0, 0.0, 0.0]
    f = 0.0
    IS = [1.0, 1.0, 1.0, 1.0]
    Ks_Zs = fill(15.475656399419266, 4)
    L = fill(0.186873685240403, 4)
    Lk = 0.017134579640331
    numn = 4
    nVG = fill(1.186873685240403, 4)
    O33 = fill(0.279610164940800, 4)
    Ohy = fill(0.096123264288629, 4)
    Osat = fill(0.459478244940800, 4)
    Pe = fill(3.912121076244333, 4)
    Qi_in = fill(0.0, 5)
    sinalp = 0.0
    Slo_pot = fill(0.0, 4)
    SN = 0.0
    T_H = fill(0.0, 4)
    T_L = fill(0.0, 4)
    V = [1.834869006521712, 1.834869006521712, 5.504607019565135, 9.174345032608560]
    Zs = [0.0, 10.0, 20.0, 50.0, 100.0]

    @testset "SPAR 1" begin
        dV = soil_moistures_rich_comp(
            V,
            Lk,
            f,
            EG,
            T_H,
            T_L,
            Qi_in,
            IS,
            1,
            Osat,
            Ohy,
            O33,
            dz,
            Ks_Zs,
            Dz,
            numn,
            L,
            Pe,
            aR,
            aT,
            alpVG,
            nVG,
            cosalp,
            sinalp,
            SN,
        )

        @test dV ≈ [-0.000048975234986, 0, 0, -0.017087673035789]
    end

    @testset "SPAR 2" begin
        dV = soil_moistures_rich_comp(
            V,
            Lk,
            f,
            EG,
            T_H,
            T_L,
            Qi_in,
            IS,
            2,
            Osat,
            Ohy,
            O33,
            dz,
            Ks_Zs,
            Dz,
            numn,
            L,
            Pe,
            aR,
            aT,
            alpVG,
            nVG,
            cosalp,
            sinalp,
            SN,
        )

        @test dV ≈ [-0.017136648270775, 0, 0, 0]
    end
end
