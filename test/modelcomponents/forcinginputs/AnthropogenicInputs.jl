using Test
using NCDatasets
using Dates
using UrbanTethysChloris.ModelComponents.ForcingInputs: initialize_anthropogenic_inputs
using ....TestUtils: load_test_netcdf

FT = Float64

hours = DateTime(2020, 1, 1, 0, 0, 0):Hour(1):DateTime(2020, 1, 1, 2, 0, 0)

filename = joinpath(mktempdir(), "test_data.nc")
ds = NCDataset(filename, "c")

# Common dimensions
defDim(ds, "hours", length(hours))
defVar(ds, "hours", hours, ("hours",))

Tbmin = 20.0
Tbmax = 25.0
defVar(ds, "Tbmin", 20.0, ())
defVar(ds, "Tbmax", 25.0, ())

Tatm = [0, 400, 296.41]

@testset "Empty netcdf" begin
    af = initialize_anthropogenic_inputs(FT, ds, Tatm)
    @test af.Tb == [Tbmin + 273.15, Tbmax + 273.15, Tatm[3]]
    @test all(af.Qf_canyon .== 0.0)
    @test all(af.Qf_roof .== 0.0)
    @test all(af.Waterf_canyonVeg .== 0.0)
    @test all(af.Waterf_canyonBare .== 0.0)
    @test all(af.Waterf_roof .== 0.0)
end

@testset "Hard-coded" begin
    Qf_canyon = [1.0, 2.0, 3.0]
    Qf_roof = [4.0, 5.0, 6.0]
    Waterf_canyonVeg = 7.0
    Waterf_canyonBare = 8.0
    Waterf_roof = [9.0, 10.0, 11.0]

    defVar(ds, "Qf_canyon", Qf_canyon, ("hours",))
    defVar(ds, "Qf_roof", Qf_roof, ("hours",))
    defVar(ds, "Waterf_canyonVeg", Waterf_canyonVeg, ())
    defVar(ds, "Waterf_canyonBare", Waterf_canyonBare, ())
    defVar(ds, "Waterf_roof", Waterf_roof, ("hours",))

    af = initialize_anthropogenic_inputs(FT, ds, Tatm)
    @test af.Tb == [Tbmin + 273.15, Tbmax + 273.15, Tatm[3]]
    @test af.Qf_canyon == Qf_canyon
    @test af.Qf_roof == Qf_roof
    @test all(af.Waterf_canyonVeg .== Waterf_canyonVeg)
    @test all(af.Waterf_canyonBare .== Waterf_canyonBare)
    @test af.Waterf_roof == Waterf_roof
end

@testset "Test file" begin
    test_ds = load_test_netcdf()
    @test_nowarn initialize_anthropogenic_inputs(FT, test_ds, Array(test_ds["Tatm"]))
end
