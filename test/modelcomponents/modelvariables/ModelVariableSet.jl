using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    ModelVariableSet, initialize_model_variable_set, TimeSlice, TimeSeries
using UrbanTethysChloris.ModelComponents.Parameters: initialize_parameter_set
using ....TestUtils: load_test_parameters

FT = Float64

input_data = load_test_parameters()

ps = initialize_parameter_set(Float64, input_data)

@testset "TimeSlice" begin
    mv_set = initialize_model_variable_set(FT, TimeSlice(), ps.soil, ps.vegetation)

    @test mv_set isa ModelVariableSet{FT,0,1}
end

@testset "TimeSeries" begin
    hours = 24
    Tatm = FT(300.0)
    AtmSpecific = FT(0.01)
    initial_value = FT(400.0)
    mv_set = initialize_model_variable_set(
        FT, TimeSeries(), Tatm, AtmSpecific, ps.soil, ps.vegetation, initial_value, hours
    )

    @test mv_set isa ModelVariableSet{FT,1,2}
end
