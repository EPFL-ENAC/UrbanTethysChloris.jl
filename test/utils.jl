using Test
using TethysChlorisCore
using UrbanTethysChloris: UrbanTethysChloris
using UrbanTethysChloris: check_extraneous_fields

FT = Float64
struct MockComponent{FT} <: AbstractModelComponent{FT}
    field1::FT
    field2::FT
    field3::FT
end

function UrbanTethysChloris.get_optional_fields(
    ::Type{MockComponent{FT}}
) where {FT<:AbstractFloat}
    return [:field1]
end

function UrbanTethysChloris.get_calculated_fields(
    ::Type{MockComponent{FT}}
) where {FT<:AbstractFloat}
    return [:field3]
end

@testset "check_extraneous_fields" begin
    # Valid case - only required fields
    data = Dict{String,Any}("field1" => 1, "field2" => 2, "field4" => 3)

    @test_nowarn check_extraneous_fields(
        MockComponent{FT}, data, ["field1", "field2", "field4"]
    )
    @test_nowarn check_extraneous_fields(
        MockComponent{FT}, data, NTuple{3,String}(("field1", "field2", "field4"))
    )
    @test_throws MethodError check_extraneous_fields(
        MockComponent{FT}, data, [:field1, :field2]
    )
    @test_throws MethodError check_extraneous_fields(
        MockComponent{FT}, data, NTuple{2,Symbol}((:field1, :field2))
    )
    @test_throws ArgumentError check_extraneous_fields(MockComponent{FT}, data, ["field1"])
    @test_throws ArgumentError check_extraneous_fields(MockComponent{FT}, data)
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

@testset "get_required_fields" begin

    # Check if the required fields are returned correctly
    required_fields = get_required_fields(MockComponent{FT})
    @test required_fields == [:field2]
end
