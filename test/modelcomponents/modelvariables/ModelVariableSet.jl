using Test
using UrbanTethysChloris.ModelComponents.ModelVariables: ModelVariableSet
using UrbanTethysChloris.ModelComponents.Parameters: initialize_soil_parameters
using ....TestUtils: load_test_parameters

FT = Float64

input_data = load_test_parameters()

soil = initialize_soil_parameters(FT, input_data["soil"])

@testset "TimeSlice" begin
    mv_set = ModelVariableSet(FT, soil)

    @test mv_set isa ModelVariableSet{FT,soil.roof.ms,soil.ground.ms}
end
