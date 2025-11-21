using Test
using UrbanTethysChloris.ModelComponents: TimeSeries
import UrbanTethysChloris.ModelComponents.Parameters: initialize_locationproperties
import UrbanTethysChloris.ModelComponents.ForcingInputs: ForcingInputSet
import ....TestUtils: load_test_parameters, load_test_netcdf

FT = Float64

@testset "ForcingInputSet" begin
    test_data = load_test_parameters()
    test_ds = load_test_netcdf()

    location_properties = initialize_locationproperties(Float64, test_data["location"])

    fis = ForcingInputSet(FT, TimeSeries(), test_ds, location_properties);

    @test fis isa ForcingInputSet{FT,1}
    @test fis[1] isa ForcingInputSet{FT,0}
end
