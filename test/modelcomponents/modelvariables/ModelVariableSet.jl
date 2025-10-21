using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    ModelVariableSet, initialize_model_variable_set
using UrbanTethysChloris.ModelComponents.Parameters: initialize_parameter_set
using ....TestUtils: load_test_parameters

FT = Float64

input_data = load_test_parameters()

ps = initialize_parameter_set(Float64, input_data)

Tatm = FT(300.0)
AtmSpecific = FT(0.01)
initial_value = FT(400.0)

@testset "ModelVariables 'scalar'" begin
    mv_set = initialize_model_variable_set(
        FT, 0, Tatm, AtmSpecific, ps.soil, ps.vegetation, initial_value
    )

    @test mv_set isa ModelVariableSet{FT,0,1}
end

@testset "ModelVariables 'model'" begin
    hours = 24
    mv_set = initialize_model_variable_set(
        FT, 1, Tatm, AtmSpecific, ps.soil, ps.vegetation, initial_value, hours
    )

    @test mv_set isa ModelVariableSet{FT,1,2}
end
