using Test
import UrbanTethysChloris.ModelComponents.Parameters: initialize_parameter_set
import ....TestUtils: load_test_parameters

input_data = load_test_parameters()

@test_nowarn initialize_parameter_set(Float64, input_data)
