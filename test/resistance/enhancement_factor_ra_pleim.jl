using Test
using UrbanTethysChloris.Resistance: enhancement_factor_ra_pleim

FT = Float64

@testset "MATLAB" begin
    disp_h = 4.902710934818196
    hPBL = 1000.0
    ra = 104.0334205197865
    Ws = 0.4
    zatm = 25.0
    zoh = 0.036417730666322
    zom = 0.364177306663220

    fconv, ra_enhanced, ra_out, LAN = enhancement_factor_ra_pleim(
        ra, zom, zoh, disp_h, zatm, Ws, hPBL
    )

    @test fconv ≈ 0.717520839172734
    @test ra_enhanced ≈ 29.387273326419383
    @test ra_out == ra
    @test LAN ≈ -3.643956285022765
end
