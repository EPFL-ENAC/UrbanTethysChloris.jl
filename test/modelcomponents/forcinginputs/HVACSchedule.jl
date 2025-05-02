using Test
using UrbanTethysChloris.ModelComponents.ForcingInputs: initialize_hvacschedule

FT = Float64

hvac_data = Dict{String,Any}(
    "Hequip" => 0.0,
    "Hpeople" => 0.0,
    "LEequip" => 0.0,
    "LEpeople" => 0.0,
    "AirConRoomFraction" => 1.0,
)

@test_nowarn initialize_hvacschedule(FT, hvac_data)
