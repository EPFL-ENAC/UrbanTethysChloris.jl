abstract type AbstractModelMetadata end
Base.@kwdef mutable struct ModelMetadata <: AbstractModelMetadata
    const ncdfpath::String
    const yamlpath::String
end

abstract type AbstractModel{FT<:AbstractFloat} end
const AM = AbstractModel

Base.@kwdef struct Model{FT<:AbstractFloat} <: AM{FT}
    variables::ModelComponents.ModelVariables.ModelVariableSet{FT}
    parameters::ModelComponents.Parameters.ParameterSet{FT}
    forcing::ModelComponents.ForcingInputs.ForcingInputSet{FT}
    metadata::AbstractModelMetadata
end

"""
    initialize_model(::Type{FT}) where {FT<:AbstractFloat}

Initialize a TethysChloris model by loading data from a NetCDF file and a YAML file.

# Arguments
- `FT`: The floating-point type to use (e.g., `Float32` or `Float64`).
- `netcdfpath`: The path to the NetCDF file containing model input data.
- `yamlpath`: The path to the YAML file containing model parameters.

# Returns
- A `Model{FT}` instance initialized with the provided data.
"""
function initialize_model(
    ::Type{FT}, netcdfpath::AbstractString, yamlpath::AbstractString
) where {FT<:AbstractFloat}
    yamldata = YAML.load_file(yamlpath; dicttype=Dict{String,Any});

    netcdfdata = NCDataset(netcdfpath);

    model = NCDataset(netcdfpath) do netcdfdata
        _initialize_model(FT, netcdfdata, yamldata, abspath(netcdfpath), abspath(yamlpath))
    end

    return model
end

function _initialize_model(
    ::Type{FT},
    netcdfdata::NCDataset,
    yamldata::Dict,
    netcdfpath::AbstractString="",
    yamlpath::AbstractString="",
) where {FT<:AbstractFloat}
    parameters = initialize_parameter_set(FT, yamldata)
    forcing = initialize_forcinginputset(FT, netcdfdata, parameters.location)
    variables = initialize_model_variable_set(
        FT,
        ModelComponents.ModelVariables.TimeSeries(),
        FT(forcing.meteorological.Tatm[1]),
        FT(forcing.meteorological.q_atm[1]),
        parameters.soil,
        parameters.vegetation,
        FT(400.0),
        length(forcing.datetime),
    )
    metadata = ModelMetadata(abspath(netcdfpath), abspath(yamlpath))

    return Model{FT}(;
        variables=variables, parameters=parameters, forcing=forcing, metadata=metadata
    )
end
