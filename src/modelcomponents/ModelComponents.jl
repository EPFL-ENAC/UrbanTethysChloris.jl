module ModelComponents

abstract type ModelDimension end
struct TimeSlice <: ModelDimension end
struct TimeSeries <: ModelDimension end
dimension_value(::TimeSlice) = 0
dimension_value(::TimeSeries) = 1
dimensionality_type(dim_value::Int) = dim_value == 0 ? TimeSlice() : TimeSeries()

export TimeSlice, TimeSeries, ModelDimension, dimension_value, dimensionality_type

include(joinpath("parameters", "Parameters.jl"))
using .Parameters

include(joinpath("forcinginputs", "ForcingInputs.jl"))
using .ForcingInputs

include(joinpath("modelvariables", "ModelVariables.jl"))
using .ModelVariables

export initialize_parameter_set, ModelVariableSet
export ParameterSet, ForcingInputSet, ModelVariableSet

end
