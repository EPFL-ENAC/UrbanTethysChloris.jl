using Test
using NCDatasets
using Dates
using UrbanTethysChloris.ModelComponents.ForcingInputs: initialize_hvacschedule
using ....TestUtils: load_test_netcdf

FT = Float64

hours = DateTime(2020, 1, 1, 0, 0, 0):Hour(1):DateTime(2020, 1, 1, 2, 0, 0)

filename = joinpath(mktempdir(), "test_data.nc")
ds = NCDataset(filename, "c")

# Common dimensions
defDim(ds, "hours", length(hours))
defVar(ds, "hours", hours, ("hours",))

@testset "Empty netcdf" begin
    hvac = initialize_hvacschedule(FT, ds)
    @test all(hvac.Hequip .== 0.0)
    @test all(hvac.Hpeople .== 0.0)
    @test all(hvac.LEequip .== 0.0)
    @test all(hvac.LEpeople .== 0.0)
    @test all(hvac.AirConRoomFraction .== 1.0)
end

@testset "Hard-coded" begin
    Hequip = [1.0, 2.0, 3.0]
    Hpeople = [4.0, 5.0, 6.0]
    LEequip = 7.0
    LEpeople = 8.0
    AirConRoomFraction = [0.1, 0.2, 0.3]

    defVar(ds, "Hequip", Hequip, ("hours",))
    defVar(ds, "Hpeople", Hpeople, ("hours",))
    defVar(ds, "LEequip", LEequip, ())
    defVar(ds, "LEpeople", LEpeople, ())
    defVar(ds, "AirConRoomFraction", AirConRoomFraction, ("hours",))

    hvac = initialize_hvacschedule(FT, ds)
    @test hvac.Hequip == Hequip
    @test hvac.Hpeople == Hpeople
    @test all(hvac.LEequip .== LEequip)
    @test all(hvac.LEpeople .== LEpeople)
    @test hvac.AirConRoomFraction == AirConRoomFraction
end

@testset "Test file" begin
    test_ds = load_test_netcdf()
    @test_nowarn initialize_hvacschedule(FT, test_ds)
end
