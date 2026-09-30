module ModelComponents

using TethysChlorisCore
using TethysChlorisCore: AbstractOutputsToSave, NoOutputs

abstract type ModelDimension end
struct TimeSlice <: ModelDimension end
struct TimeSeries <: ModelDimension end
dimension_value(::TimeSlice) = 0
dimension_value(::TimeSeries) = 1
dimensionality_type(dim_value::Int) = dim_value == 0 ? TimeSlice() : TimeSeries()

export TimeSlice, TimeSeries, ModelDimension, dimension_value, dimensionality_type

# In increasing order, Plot outputs are mandatory for plotting, the others are supersets
# of the previous ones, providing outputs with an increasing level of detail.
struct PlotOutputs <: AbstractOutputsToSave end
struct EssentialOutputs <: AbstractOutputsToSave end
struct ExtendedEnergyClimateOutputs <: AbstractOutputsToSave end
struct ExtendedOutputs <: AbstractOutputsToSave end

const no_outputs = NoOutputs()
const plot_outputs = PlotOutputs()
const essential_outputs = EssentialOutputs()
const extended_energy_climate_outputs = ExtendedEnergyClimateOutputs()
const extended_outputs = ExtendedOutputs()

export no_outputs,
    plot_outputs, essential_outputs, extended_energy_climate_outputs, extended_outputs

TethysChlorisCore.decrease(::Type{PlotOutputs}) = NoOutputs
TethysChlorisCore.decrease(::Type{EssentialOutputs}) = PlotOutputs
TethysChlorisCore.decrease(::Type{ExtendedEnergyClimateOutputs}) = EssentialOutputs
TethysChlorisCore.decrease(::Type{ExtendedOutputs}) = ExtendedEnergyClimateOutputs

export PlotOutputs, EssentialOutputs, ExtendedEnergyClimateOutputs, ExtendedOutputs

include(joinpath("parameters", "Parameters.jl"))
using .Parameters

include(joinpath("forcinginputs", "ForcingInputs.jl"))
using .ForcingInputs

include(joinpath("modelvariables", "ModelVariables.jl"))
using .ModelVariables

export initialize_parameter_set, ModelVariableSet
export ParameterSet, ForcingInputSet, ModelVariableSet
export update!

end
