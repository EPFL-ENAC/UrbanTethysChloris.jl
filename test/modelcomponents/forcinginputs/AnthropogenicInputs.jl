using Test
using NCDatasets
using Dates
using UrbanTethysChloris.ModelComponents: TimeSeries
using UrbanTethysChloris.ModelComponents.ForcingInputs: AnthropogenicInputs
using ....TestUtils: load_test_netcdf

FT = Float64

Tbmin = FT(20.0)
Tbmax = FT(25.0)
Tatm = FT.([0, 400, 296.41])
hours = length(Tatm)

@testset "AnthropogenicInputs" begin
    @testset "Empty netcdf" begin
        af = mktempdir() do tempdir
            filename = joinpath(tempdir, "empty.nc")
            ds = NCDataset(filename, "c")

            defDim(ds, "hours", hours)
            defVar(ds, "Tbmin", Tbmin, ())
            defVar(ds, "Tbmax", Tbmax, ())

            AnthropogenicInputs(FT, TimeSeries(), ds, Tatm)
        end

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

        af = mktempdir() do tempdir
            filename = joinpath(tempdir, "hardcoded.nc")

            ds = NCDataset(filename, "c")

            defDim(ds, "hours", hours)
            defVar(ds, "Tbmin", Tbmin, ())
            defVar(ds, "Tbmax", Tbmax, ())

            defVar(ds, "Qf_canyon", Qf_canyon, ("hours",))
            defVar(ds, "Qf_roof", Qf_roof, ("hours",))
            defVar(ds, "Waterf_canyonVeg", Waterf_canyonVeg, ())
            defVar(ds, "Waterf_canyonBare", Waterf_canyonBare, ())
            defVar(ds, "Waterf_roof", Waterf_roof, ("hours",))

            AnthropogenicInputs(FT, TimeSeries(), ds, Tatm)
        end

        @test af.Tb == [Tbmin + 273.15, Tbmax + 273.15, Tatm[3]]
        @test af.Qf_canyon == Qf_canyon
        @test af.Qf_roof == Qf_roof
        @test all(af.Waterf_canyonVeg .== Waterf_canyonVeg)
        @test all(af.Waterf_canyonBare .== Waterf_canyonBare)
        @test af.Waterf_roof == Waterf_roof
    end

    @testset "Test file" begin
        ds = load_test_netcdf();

        ai = AnthropogenicInputs(FT, TimeSeries(), ds, Vector(ds["Tatm"]));
        hours = ds.dim["hours"]

        @test ai isa AnthropogenicInputs{FT,1}

        # Test initial values
        @test ai.Qf_canyon[1] == ds["Qf_canyon"][1]

        # Test all fields are accessible and have correct dimensions
        for field in fieldnames(AnthropogenicInputs)
            @test isa(getproperty(ai, field), Array{FT,1})
            @test size(getproperty(ai, field)) == (hours,)
        end

        @test ai[1] isa AnthropogenicInputs{FT,0}
    end
end
