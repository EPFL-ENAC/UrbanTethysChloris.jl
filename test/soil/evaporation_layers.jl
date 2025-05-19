using Test
using UrbanTethysChloris.Soil: evaporation_layers

@testset "MATLAB" begin
    # Test case 1: Basic functionality with valid inputs
    Zs = [0.0, 100.0, 200.0, 300.0]
    Zdes = 150.0
    result = evaporation_layers(Zs, Zdes)
    @test result ≈ [2/3, 1/3, 0.0]

    # Test case 2: Desorption depth equals layer boundary
    Zs = [0.0, 50.0, 100.0]
    Zdes = 50.0
    result = evaporation_layers(Zs, Zdes)
    @test result ≈ [1.0, 0.0]

    # Test case 3: Desorption depth too shallow (less than first layer)
    Zs = [10.0, 20.0, 30.0]
    Zdes = 5.0
    result = evaporation_layers(Zs, Zdes)
    @test all(result .== 0.0)

    # Test case 4: Desorption depth too deep (more than last layer)
    Zs = [0.0, 100.0, 200.0]
    Zdes = 250.0
    result = evaporation_layers(Zs, Zdes)
    @test all(result .== 0.0)

    # Test case 5: Edge case with minimum valid layers
    Zs = [0.0, 100.0]
    Zdes = 50.0
    result = evaporation_layers(Zs, Zdes)
    @test length(result) == 1
    @test result[1] == 1.0
end
