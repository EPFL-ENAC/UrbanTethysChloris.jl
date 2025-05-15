using Test
using NCDatasets
using Dates
using UrbanTethysChloris.ModelComponents.ForcingInputs:
    SunPositionInputs, initialize_sunposition_inputs
using UrbanTethysChloris.ModelComponents.Parameters: initialize_locationproperties
using ....TestUtils: load_test_netcdf, load_test_parameters

FT = Float64

@testset "Hard-coded" begin
    hours = DateTime(2020, 1, 1, 0, 0, 0):Hour(1):DateTime(2020, 1, 1, 2, 0, 0)

    filename = joinpath(mktempdir(), "test_data.nc")
    ds = NCDataset(filename, "c")

    # Common dimensions
    defDim(ds, "hours", length(hours))
    defVar(ds, "hours", hours, ("hours",))

    defVar(ds, "t_bef", 0.5, ())
    defVar(ds, "t_aft", 0.5, ())

    location_data = Dict{String,Any}(
        "phi" => 25.3, "lambda" => 55.4, "theta_canyon" => deg2rad(45.0), "DeltaGMT" => 4.0
    )
    location = initialize_locationproperties(FT, location_data)

    @test_nowarn initialize_sunposition_inputs(FT, ds, Array(ds["hours"]), location)
    close(ds)
end

@testset "Test file" begin
    yaml_data = load_test_parameters()
    location = initialize_locationproperties(FT, yaml_data["location"])
    test_ds = load_test_netcdf()
    @test_nowarn initialize_sunposition_inputs(
        FT, test_ds, Array(test_ds["datetime"]), location
    )
end
