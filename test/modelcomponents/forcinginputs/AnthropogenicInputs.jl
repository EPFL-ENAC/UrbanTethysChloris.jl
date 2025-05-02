using Test
using NCDatasets
using Dates
using UrbanTethysChloris.ModelComponents.ForcingInputs:
    AnthropogenicInputs, initialize_anthropogenic_inputs
using ....TestUtils: load_test_netcdf

FT = Float64

@testset "Hard-coded" begin
    hours = DateTime(2020, 1, 1, 0, 0, 0):Hour(1):DateTime(2020, 1, 1, 2, 0, 0)

    filename = joinpath(mktempdir(), "test_data.nc")
    ds = NCDataset(filename, "c")

    # Common dimensions
    defDim(ds, "hours", length(hours))
    defVar(ds, "hours", hours, ("hours",))

    defVar(ds, "Tbmin", 20.0, ())
    defVar(ds, "Tbmax", 25.0, ())

    @test_nowarn initialize_anthropogenic_inputs(FT, ds, [0, 400, 296.41])
end

@testset "Test file" begin
    test_ds = load_test_netcdf()
    @test_nowarn initialize_anthropogenic_inputs(FT, test_ds, Array(test_ds["Tatm"]))
end
