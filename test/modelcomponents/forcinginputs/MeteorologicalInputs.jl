using Test
using UrbanTethysChloris.ModelComponents: TimeSeries
using UrbanTethysChloris.ModelComponents.ForcingInputs: MeteorologicalInputs
using Dates
using ....TestUtils: load_test_netcdf

FT = Float64

@testset "MeteorologicalInputs" begin
    ds = load_test_netcdf();
    theta_Z = zeros(FT, ds.dim["hours"])
    hours = ds.dim["hours"]

    mi = MeteorologicalInputs(FT, TimeSeries(), ds, theta_Z);

    @test mi isa MeteorologicalInputs{FT,1}

    # Test initial values
    @test mi.LWR_in[1] == ds["LWR_in"][1]

    scalar_meteo_fields = [:Zatm, :Catm_CO2, :Catm_O2, :SunDSM_MRT]

    # Test all fields are accessible and have correct dimensions
    for field in fieldnames(MeteorologicalInputs)
        if field in scalar_meteo_fields
            @test isa(getproperty(mi, field), FT)
        elseif field == :datetime
            @test isa(getproperty(mi, field), Array{DateTime,1})
            @test size(getproperty(mi, field)) == (hours,)
        else
            @test isa(getproperty(mi, field), Array{FT,1})
            @test size(getproperty(mi, field)) == (hours,)
        end
    end

    @test mi[1] isa MeteorologicalInputs{FT,0}
end
