using Test
using NCDatasets
using Dates
using UrbanTethysChloris.ModelComponents: TimeSeries
using UrbanTethysChloris.ModelComponents.ForcingInputs: SunPositionInputs
using UrbanTethysChloris.ModelComponents.Parameters: initialize_locationproperties
using ....TestUtils: load_test_netcdf, load_test_parameters

FT = Float64

@testset "SunPositionInputs" begin
    @testset "Hard-coded" begin
        hours = DateTime(2020, 1, 1, 0, 0, 0):Hour(1):DateTime(2020, 1, 1, 2, 0, 0)

        sp = mktempdir() do tempdir
            filename = joinpath(tempdir, "hardcoded.nc")
            ds = NCDataset(filename, "c")

            defDim(ds, "hours", length(hours))
            defVar(ds, "hours", hours, ("hours",))

            defVar(ds, "t_bef", 0.5, ())
            defVar(ds, "t_aft", 0.5, ())

            location_data = Dict{String,Any}(
                "phi" => 25.3,
                "lambda" => 55.4,
                "theta_canyon" => deg2rad(45.0),
                "DeltaGMT" => 4.0,
            )
            location = initialize_locationproperties(FT, location_data)

            SunPositionInputs(FT, TimeSeries(), ds, Array(ds["hours"]), location)
        end

        @test sp isa SunPositionInputs{FT,1}
        @test sp.t_bef == 0.5
        @test sp.t_aft == 0.5
    end

    @testset "Test file" begin
        yaml_data = load_test_parameters()
        location = initialize_locationproperties(FT, yaml_data["location"])

        ds = load_test_netcdf();
        hours = ds.dim["hours"]

        sp = SunPositionInputs(FT, TimeSeries(), ds, Array(ds["datetime"]), location);

        @test sp isa SunPositionInputs{FT,1}

        scalar_sunpos_fields = [:t_bef, :t_aft]

        # Test all fields are accessible and have correct dimensions
        for field in fieldnames(SunPositionInputs)
            if field in scalar_sunpos_fields
                @test isa(getproperty(sp, field), FT)
            else
                @test isa(getproperty(sp, field), Array{FT,1})
                @test size(getproperty(sp, field)) == (hours,)
            end
        end

        @test sp[1] isa SunPositionInputs{FT,0}
    end
end
