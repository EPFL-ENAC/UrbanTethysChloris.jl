using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    BuildingEnergyModelVariables, initialize_building_energy_model_variables

FT = Float64

@testset "BuildingEnergyModelVariables scalar (N=0)" begin
    bem_vars = initialize_building_energy_model_variables(FT, 0)

    # Test structure
    @test bem_vars isa BuildingEnergyModelVariables{FT,0}

    # Test fields initialized to zero
    @test bem_vars.HBinRoof == 0.0
    @test bem_vars.SWRabsGround == 0.0
    @test bem_vars.AC_on == 0.0

    # Test all fields are accessible
    for field in fieldnames(BuildingEnergyModelVariables)
        @test isa(getproperty(bem_vars, field), FT)
    end
end

@testset "BuildingEnergyModelVariables vector (N=1)" begin
    hours = 24
    tatm = FT(300.0)  # Initial temperature
    qatm = FT(0.01)   # Initial humidity

    bem_vars = initialize_building_energy_model_variables(FT, 1, hours, tatm, qatm)

    # Test structure
    @test bem_vars isa BuildingEnergyModelVariables{FT,1}

    # Test field dimensions
    @test size(bem_vars.Tbin) == (hours,)
    @test size(bem_vars.qbin) == (hours,)

    # Test initial values
    @test bem_vars.Tbin[1] == tatm
    @test bem_vars.qbin[1] == qatm

    # Test all fields are accessible and have correct dimensions
    for field in fieldnames(BuildingEnergyModelVariables)
        @test isa(getproperty(bem_vars, field), Array{FT,1})
        @test size(getproperty(bem_vars, field)) == (hours,)
    end
end
