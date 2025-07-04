using SafeTestsets

@safetestset "AirHumidity2m" begin
    include("air_humidity_2m.jl")
end

@safetestset "AirHumidity2mOutput" begin
    include("air_humidity_2m_output.jl")
end

@safetestset "CalculateT2m" begin
    include("calculate_t2m.jl")
end

@safetestset "HeatFluxCanyon" begin
    include("heat_flux_canyon.jl")
end

@safetestset "HeatFlux_ground" begin
    include("heat_flux_ground.jl")
end

@safetestset "HeatFluxRoof" begin
    include("heat_flux_roof.jl")
end

@safetestset "HeatFluxWall" begin
    include("heat_flux_wall.jl")
end
