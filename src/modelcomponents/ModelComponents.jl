module ModelComponents

using Infiltrator

abstract type ModelDimension end
struct TimeSlice <: ModelDimension end
struct TimeSeries <: ModelDimension end
dimension_value(::TimeSlice) = 0
dimension_value(::TimeSeries) = 1
dimensionality_type(dim_value::Int) = dim_value == 0 ? TimeSlice() : TimeSeries()

export TimeSlice, TimeSeries, ModelDimension, dimension_value, dimensionality_type

abstract type AbstractOutputsToSave end

# In increasing order, Plot outputs are mandatory for plotting, the others are supersets
# of the previous ones, providing outputs with an increasing level of detail.
struct NoOutputs <: AbstractOutputsToSave end
struct PlotOutputs <: AbstractOutputsToSave end
struct EssentialOutputs <: AbstractOutputsToSave end
struct ExtendedEnergyClimateOutputs <: AbstractOutputsToSave end
struct ExtendedOutputs <: AbstractOutputsToSave end

const no_outputs = NoOutputs()
const plot_outputs = PlotOutputs()
const essential_outputs = EssentialOutputs()
const extended_energy_climate_outputs = ExtendedEnergyClimateOutputs()
const extended_outputs = ExtendedOutputs()

decrease(::Type{PlotOutputs}) = NoOutputs
decrease(::Type{EssentialOutputs}) = PlotOutputs
decrease(::Type{ExtendedEnergyClimateOutputs}) = EssentialOutputs
decrease(::Type{ExtendedOutputs}) = ExtendedEnergyClimateOutputs

function outputs_to_save(::Type{T}) where {T}
    return fieldnames(T)
end

"""
    outputs_to_save(::Type{T}, ::Type{O}) where {T, O <: AbstractOutputsToSave}

Get the outputs to save for type `T` corresponding to the outputs to save defined by type
`AbstractOutputsToSave`. By default, returns an empty tuple. This function is expected to be
overloaded for specific types.

# Arguments
- `T::Type{T}`: The model component type
- `::Type{AbstractOutputsToSave}`: The outputs to save type

# Returns
- `Tuple`: Tuple of output field names
"""
function outputs_to_save(::Type{T}, ::Type{O}) where {T,O<:AbstractOutputsToSave}
    return ()
end

"""
    accessors(T::Type{T}, O::Type{O}) where {T,O}

Create a nested dictionary of accessor functions for type `T` corresponding to the outputs
to save defined by type `O`. The function is recursive, building the dictionary from the
bottom up by decreasing the output type until reaching `NoOutputs`.

# Arguments
- `T::Type{T}`: The model component type
- `O::Type{O}`: The outputs to save type

# Returns
- `Dict{Symbol,Dict{Symbol,Function}}`: Nested dictionary of accessor functions
"""
function accessors(::Type{T}, ::Type{O}) where {T,O}
    fns = ModelComponents.outputs_to_save(T, O)

    base = accessors(T, decrease(O))
    if !isempty(fns)
        merge!(base, create_nested_accessor_dict(T, fns))
    end
    return base
end

function accessors(::Type{T}, ::Type{NoOutputs}) where {T}
    return Dict{Symbol,Dict{Symbol,Function}}()
end

"""
    create_accessor_dict(::Type{T}) where {T}

Create a dictionary of accessor functions for type `T`.
"""
function create_accessor_dict(::Type{T}) where {T}
    return Dict{Symbol,Function}(fn => (x -> getfield(x, fn)) for fn in fieldnames(T))
end

"""
    create_nested_accessor_dict(::Type{T}, fns::NTuple=fieldnames(T)) where {T}
Create a nested dictionary of accessor functions for type `T` for the specified fields `fns`.

# Arguments
- `T::Type{T}`: The model component type
- `fns::NTuple`: The fields to create accessors for (default: all fields of `T`)

# Returns
- `Dict{Symbol,Dict{Symbol,Function}}`: Nested dictionary of accessor functions
"""
function create_nested_accessor_dict(::Type{T}, fns::NTuple=fieldnames(T)) where {T}
    return Dict{Symbol,Dict{Symbol,Function}}(
        fn => Dict{Symbol,Function}(
            sub_fn => (x -> getfield(getfield(x, fn), sub_fn)) for
            sub_fn in fieldnames(fieldtype(T, fn))
        ) for fn in fns
    )
end

export outputs_to_save,
    EssentialOutputs, PlotOutputs, ExtendedEnergyClimateOutputs, ExtendedOutputs

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
