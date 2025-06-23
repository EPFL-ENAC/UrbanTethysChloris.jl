using SafeTestsets

@safetestset "ACHeatingModule" begin
    include("ac_heating_module.jl")
end

@safetestset "ACHeatingTurnOnOff" begin
    include("ac_heating_turn_on_off.jl")
end
