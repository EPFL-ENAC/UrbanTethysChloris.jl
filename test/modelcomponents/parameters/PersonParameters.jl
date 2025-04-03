using Test
using UrbanTethysChloris.ModelComponents.Parameters:
    PersonParameters, initialize_person_parameters

FT = Float64

person_params = Dict{String,Any}(
    "PositionPx" => 5.78 / 2,
    "PositionPz" => 1.1,
    "PersonWidth" => 0.06 / 2,
    "PersonHeight" => 0.22 / 2,
    "HeightWind" => 1.1,
)

@testset "PersonParameters initialization" begin
    person = initialize_person_parameters(FT, person_params)

    @test person.PositionPx == FT(5.78 / 2)
    @test person.PositionPz == FT(1.1)
    @test person.PersonWidth == FT(0.06 / 2)
    @test person.PersonHeight == FT(0.22 / 2)
    @test person.HeightWind == FT(1.1)
end
