using Test
using UrbanTethysChloris.Resistance: backcalculate_obhukov_length

FT = Float64

@testset "MATLAB" begin
    disp_h = 4.902710934818196
    ra = 104.0334205197865
    u = 0.4
    zatm = 25.0
    zoh = 0.036417730666322
    zom = 0.364177306663220

    # Call the function
    LAN, Diff_ra = backcalculate_obhukov_length(ra, zom, zoh, disp_h, zatm, u)

    @test LAN ≈ -3.643956285022765
end
