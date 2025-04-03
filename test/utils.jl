using Test
using TethysChlorisCore
using UrbanTethysChloris: check_extraneous_fields

FT = Float64
struct MockComponent{FT} <: AbstractModelComponent{FT} end

# Define required fields for mock component
function TethysChlorisCore.get_required_fields(
    ::Type{MockComponent{FT}}
) where {FT<:AbstractFloat}
    return [:field1, :field2]
end

@testset "check_extraneous_fields" begin
    # Valid case - only required fields
    data = Dict{String,Any}("field1" => 1, "field2" => 2)

    @test_nowarn check_extraneous_fields(MockComponent{FT}, data, ["field1", "field2"])
    @test_nowarn check_extraneous_fields(
        MockComponent{FT}, data, NTuple{2,String}(("field1", "field2"))
    )
    @test_throws MethodError check_extraneous_fields(
        MockComponent{FT}, data, [:field1, :field2]
    )
    @test_throws MethodError check_extraneous_fields(
        MockComponent{FT}, data, NTuple{2,Symbol}((:field1, :field2))
    )
    @test_throws ArgumentError check_extraneous_fields(MockComponent{FT}, data, ["field1"])
end

@testset "validate_fields" begin
    # Valid case
    valid_data = Dict{String,Any}("field1" => 1, "field2" => 2)
    @test_nowarn TethysChlorisCore.validate_fields(MockComponent{FT}, valid_data)

    # Invalid case
    invalid_data = Dict{String,Any}("field1" => 1, "field2" => 2, "extra_field" => 3)
    @test_throws ArgumentError TethysChlorisCore.validate_fields(
        MockComponent{FT}, invalid_data
    )
end
