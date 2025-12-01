using Test
using TethysChlorisCore
using UrbanTethysChloris.ModelComponents.ModelVariables:
    TempVec, TempDamp, MRT, ThermalComfort, TemperatureVariables

FT = Float64

@testset "Subsets" begin
    @testset "TempVec" begin
        tempvec = TempVec(FT)
        @test tempvec isa TempVec{FT}

        for field in fieldnames(TempVec)
            @test isa(getproperty(tempvec, field), FT)
            @test getproperty(tempvec, field) == 0
        end
    end

    @testset "TempDamp" begin
        tempdamp = TempDamp(FT)
        @test tempdamp isa TempDamp{FT}

        for field in fieldnames(TempDamp)
            @test isa(getproperty(tempdamp, field), FT)
            @test getproperty(tempdamp, field) == 0
        end
    end

    @testset "MRT" begin
        mrt = MRT(FT)
        @test mrt isa MRT{FT}

        for field in fieldnames(MRT)
            @test isa(getproperty(mrt, field), FT)
            @test getproperty(mrt, field) == 0
        end
    end

    @testset "ThermalComfort" begin
        tc = ThermalComfort(FT)
        @test tc isa ThermalComfort{FT}

        for field in fieldnames(ThermalComfort)
            @test isa(getproperty(tc, field), FT)
            @test getproperty(tc, field) == 0
        end
    end
end

@testset "TemperatureVariables" begin
    temp_vars = TemperatureVariables(FT)

    # Test structure
    @test temp_vars isa TemperatureVariables{FT}

    # Test components
    @test temp_vars.tempvec isa TempVec{FT}
    @test temp_vars.tempdamp isa TempDamp{FT}
    @test temp_vars.mrt isa MRT{FT}
    @test temp_vars.thermalcomfort isa ThermalComfort{FT}

    # Test field access for scalar case
    @test temp_vars.tempvec.TRoofImp === 0.0
    @test temp_vars.thermalcomfort.UTCI === 0.0
end
