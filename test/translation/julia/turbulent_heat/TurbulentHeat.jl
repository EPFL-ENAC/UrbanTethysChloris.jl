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
