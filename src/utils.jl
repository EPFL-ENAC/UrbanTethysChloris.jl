"""
    check_extraneous_fields(::Type{T}, data::Dict{String,Any}) where {T<:AbstractModelComponent}

Check if the input dictionary contains any keys that are not part of the required fields for type T.

# Arguments
- `T`: Type of model component to check fields against
- `data`: Dictionary containing fields to validate
- `required_fields`: List of required fields for the model component. This is typically obtained from `get_required_fields(T)` or `fieldnames(T)`

# Throws
- `ArgumentError`: If any extraneous keys are found in the data dictionary
"""
function check_extraneous_fields(
    ::Type{T},
    data::Dict{String,Any},
    required_fields::Union{Vector{String},NTuple{N,String}},
) where {T<:AbstractModelComponent,N}
    for key in keys(data)
        if key ∉ required_fields
            throw(ArgumentError("Extraneous key: $key"))
        end
    end
end

"""
    TethysChlorisCore.validate_fields(::Type{T}, data::Dict{String,Any}) where {T<:AbstractModelComponent}

Validate that a dictionary contains only the required fields for a given model component type.

# Arguments
- `T`: Type of model component to validate fields against
- `data`: Dictionary containing fields to validate

# Throws
- `ArgumentError`: If any extraneous fields are found in the data dictionary
"""
function TethysChlorisCore.validate_fields(
    ::Type{T}, data::Dict{String,Any}
) where {T<:AbstractModelComponent}
    return check_extraneous_fields(T, data, String.(get_required_fields(T)))
end
