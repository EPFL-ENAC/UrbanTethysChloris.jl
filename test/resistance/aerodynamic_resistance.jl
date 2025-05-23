using Test
using UrbanTethysChloris.Resistance: aerodynamic_resistance

FT = Float64

@testset "MATLAB" begin
    disp_h = 0.1
    ea = 763.6587960785345
    es = 807.2503129794234
    Pre = 948.7232105699336
    Ta = 3.899999999999977
    Ts = 3.899999999999977
    Ws = 0.4
    zatm = 18.5
    zoh = 0.001845
    zom = 0.01845

    ra = aerodynamic_resistance(Ta, Ts, Pre, zatm, disp_h, zom, zoh, Ws, ea, es)

    @test ra ≈ 725.1732619888475
end
