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

function check_extraneous_fields(
    ::Type{T}, data::Dict{String,Any}
) where {T<:AbstractModelComponent}
    return check_extraneous_fields(T, data, String.(fieldnames(T)))
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
    return check_extraneous_fields(T, data)
end

function TethysChlorisCore.get_required_fields(::Type{T}) where {T<:AbstractModelComponent}
    all_fields = Set(fieldnames(T))
    optional_fields = Set(get_optional_fields(T))
    calculated_fields = Set(get_calculated_fields(T))
    return collect(setdiff(setdiff(all_fields, optional_fields), calculated_fields))
end

"""
    get_optional_fields(::Type{T}) where {T<:AbstractModelComponent}

Get a list of optional fields for a given model component type. Optional fields are those that can be omitted when creating
a model component. By default, returns an empty list. Components should override this method if they have optional fields.

# Arguments
- `T`: Type of model component to get optional fields for

# Returns
- `Vector{Symbol}`: List of optional field names
"""
function get_optional_fields(::Type{T}) where {T<:AbstractModelComponent}
    return Symbol[]
end

"""
    get_calculated_fields(::Type{T}) where {T<:AbstractModelComponent}

Get a list of calculated fields for a given model component type. Calculated fields are those that are computed based on
other fields and should not be provided when creating a model component. By default, returns an empty list. Components
should override this method if they have calculated fields.

# Arguments
- `T`: Type of model component to get calculated fields for

# Returns
- `Vector{Symbol}`: List of calculated field names
"""
function get_calculated_fields(::Type{T}) where {T<:AbstractModelComponent}
    return Symbol[]
end
