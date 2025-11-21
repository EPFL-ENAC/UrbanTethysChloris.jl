using Test
using NCDatasets
using Dates
using UrbanTethysChloris.ModelComponents: TimeSeries
using UrbanTethysChloris.ModelComponents.ForcingInputs: HVACSchedule
using ....TestUtils: load_test_netcdf

FT = Float64

hours = 3

@testset "HVACSchedule" begin
    @testset "Empty netcdf" begin
        hvac = mktempdir() do tempdir
            filename = joinpath(tempdir, "empty.nc")
            ds = NCDataset(filename, "c")

            defDim(ds, "hours", hours)

            HVACSchedule(FT, TimeSeries(), ds)
        end

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

        hvac = mktempdir() do tempdir
            filename = joinpath(tempdir, "hardcoded.nc")
            ds = NCDataset(filename, "c")

            defDim(ds, "hours", hours)

            defVar(ds, "Hequip", Hequip, ("hours",))
            defVar(ds, "Hpeople", Hpeople, ("hours",))
            defVar(ds, "LEequip", LEequip, ())
            defVar(ds, "LEpeople", LEpeople, ())
            defVar(ds, "AirConRoomFraction", AirConRoomFraction, ("hours",))

            HVACSchedule(FT, TimeSeries(), ds)
        end

        @test hvac.Hequip == Hequip
        @test hvac.Hpeople == Hpeople
        @test all(hvac.LEequip .== LEequip)
        @test all(hvac.LEpeople .== LEpeople)
        @test hvac.AirConRoomFraction == AirConRoomFraction
    end

    @testset "Test file" begin
        ds = load_test_netcdf();

        hvac = HVACSchedule(FT, TimeSeries(), ds);
        hours = ds.dim["hours"]

        @test hvac isa HVACSchedule{FT,1}

        # Test all fields are accessible and have correct dimensions
        for field in fieldnames(HVACSchedule)
            @test isa(getproperty(hvac, field), Array{FT,1})
            @test size(getproperty(hvac, field)) == (hours,)
        end

        @test hvac[1] isa HVACSchedule{FT,0}
    end
end
