using Test
import UrbanTethysChloris.ModelComponents.Parameters: initialize_locationproperties
import UrbanTethysChloris.ModelComponents.ForcingInputs: initialize_forcinginputset
import ....TestUtils: load_test_parameters, load_test_netcdf

test_data = load_test_parameters()
test_ds = load_test_netcdf()

location_properties = initialize_locationproperties(Float64, test_data["location"])

@test_nowarn initialize_forcinginputset(Float64, test_ds, location_properties)
