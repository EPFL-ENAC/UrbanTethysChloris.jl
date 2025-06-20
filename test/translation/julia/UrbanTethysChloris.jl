using SafeTestsets

@safetestset "resistance" begin
    include(joinpath("resistance", "Resistance.jl"))
end

@safetestset "soil" begin
    include(joinpath("soil", "Soil.jl"))
end

@safetestset "turbulent_heat" begin
    include(joinpath("turbulent_heat", "TurbulentHeat.jl"))
end

@safetestset "water" begin
    include(joinpath("water", "Water.jl"))
end
