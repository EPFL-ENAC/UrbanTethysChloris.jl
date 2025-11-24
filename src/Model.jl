abstract type AbstractModelMetadata end
Base.@kwdef mutable struct ModelMetadata <: AbstractModelMetadata
    const ncdfpath::String
    const yamlpath::String
end

abstract type AbstractModel{FT<:AbstractFloat} end
const AM = AbstractModel

abstract type Abstract1PModel{FT<:AbstractFloat} <: AbstractModel{FT} end
Base.@kwdef mutable struct Model{FT<:AbstractFloat} <: AM{FT}
    variables::ModelComponents.ModelVariables.AbstractModelVariableSet{FT}
    const parameters::ModelComponents.Parameters.ParameterSet{FT}
    forcing::ModelComponents.ForcingInputs.ForcingInputSet{FT,0}
    const metadata::AbstractModelMetadata
end

"""
    Model(::Type{FT}) where {FT<:AbstractFloat}

Initialize a TethysChloris model by loading data from a NetCDF file and a YAML file.

# Arguments
- `FT`: The floating-point type to use (e.g., `Float32` or `Float64`).
- `netcdfpath`: The path to the NetCDF file containing model input data.
- `yamlpath`: The path to the YAML file containing model parameters.

# Returns
- A `Model{FT, MR, MG}` instance initialized with the provided data.
"""
function create_model(
    ::Type{FT}, netcdfpath::AbstractString, yamlpath::AbstractString
) where {FT<:AbstractFloat}
    yamldata = YAML.load_file(yamlpath; dicttype=Dict{String,Any});

    model, forcing = NCDataset(netcdfpath) do netcdfdata
        _create_model(FT, netcdfdata, yamldata, abspath(netcdfpath), abspath(yamlpath))
    end

    return model, forcing
end

function _create_model(
    ::Type{FT},
    netcdfdata::NCDataset,
    yamldata::Dict,
    netcdfpath::AbstractString="",
    yamlpath::AbstractString="",
) where {FT<:AbstractFloat}
    parameters = initialize_parameter_set(FT, yamldata)
    forcing = ForcingInputSet(FT, TimeSeries(), netcdfdata, parameters.location)
    variables = ModelVariableSet(FT, parameters.soil)
    metadata = ModelMetadata(abspath(netcdfpath), abspath(yamlpath))

    model = Model{FT}(;
        variables=variables, parameters=parameters, forcing=forcing[1], metadata=metadata
    )

    return model, forcing
end
