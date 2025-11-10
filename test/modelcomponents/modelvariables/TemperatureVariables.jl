using Test
using TethysChlorisCore
using UrbanTethysChloris.ModelComponents.ModelVariables:
    TempVec,
    TempDamp,
    MRT,
    ThermalComfort,
    initialize_tempvec,
    initialize_tempdamp,
    initialize_mrt,
    initialize_thermal_comfort,
    TemperatureVariables,
    initialize_temperature_variables,
    TimeSlice,
    TimeSeries

FT = Float64

hours = 24

@testset "TempVec" begin
    @testset "TimeSlice" begin
        tempvec = initialize_tempvec(FT, TimeSlice())
        @test tempvec isa TempVec{FT,0}

        for field in fieldnames(TempVec)
            @test isa(getproperty(tempvec, field), FT)
            @test getproperty(tempvec, field) == 0
        end
    end

    @testset "TimeSeries" begin
        tempvec = initialize_tempvec(FT, TimeSeries(), hours)
        @test tempvec isa TempVec{FT,1}

        for field in fieldnames(TempVec)
            @test isa(getproperty(tempvec, field), Array{FT,1})
            @test size(getproperty(tempvec, field)) == (hours,)
            @test getproperty(tempvec, field)[1] == 0
        end
    end
end

@testset "TempDamp" begin
    @testset "TimeSlice" begin
        tempdamp = initialize_tempdamp(FT, TimeSlice())
        @test tempdamp isa TempDamp{FT,0}

        for field in fieldnames(TempDamp)
            @test isa(getproperty(tempdamp, field), FT)
            @test getproperty(tempdamp, field) == 0
        end
    end

    @testset "TimeSeries" begin
        tempdamp = initialize_tempdamp(FT, TimeSeries(), hours)
        @test tempdamp isa TempDamp{FT,1}

        for field in fieldnames(TempDamp)
            @test isa(getproperty(tempdamp, field), Array{FT,1})
            @test size(getproperty(tempdamp, field)) == (hours,)
            @test getproperty(tempdamp, field)[1] == 0
        end
    end
end

@testset "MRT" begin
    @testset "TimeSlice" begin
        mrt = initialize_mrt(FT, TimeSlice())
        @test mrt isa MRT{FT,0}

        for field in fieldnames(MRT)
            @test isa(getproperty(mrt, field), FT)
            @test getproperty(mrt, field) == 0
        end
    end

    @testset "TimeSeries" begin
        mrt = initialize_mrt(FT, TimeSeries(), hours)
        @test mrt isa MRT{FT,1}

        for field in fieldnames(MRT)
            @test isa(getproperty(mrt, field), Array{FT,1})
            @test size(getproperty(mrt, field)) == (hours,)
            @test getproperty(mrt, field)[1] == 0
        end
    end
end

@testset "ThermalComfort" begin
    @testset "TimeSlice" begin
        tc = initialize_thermal_comfort(FT, TimeSlice())
        @test tc isa ThermalComfort{FT,0}

        for field in fieldnames(ThermalComfort)
            @test isa(getproperty(tc, field), FT)
            @test getproperty(tc, field) == 0
        end
    end

    @testset "TimeSeries" begin
        tc = initialize_thermal_comfort(FT, TimeSeries(), hours)
        @test tc isa ThermalComfort{FT,1}

        for field in fieldnames(ThermalComfort)
            @test isa(getproperty(tc, field), Array{FT,1})
            @test size(getproperty(tc, field)) == (hours,)
            @test getproperty(tc, field)[1] == 0
        end
    end
end

@testset "TemperatureVariables" begin
    @testset "TimeSlice" begin
        temp_vars = initialize_temperature_variables(FT, TimeSlice())

        # Test structure
        @test temp_vars isa TemperatureVariables{FT,0}

        # Test components
        @test temp_vars.tempvec isa TempVec{FT,0}
        @test temp_vars.tempdamp isa TempDamp{FT,0}
        @test temp_vars.mrt isa MRT{FT,0}
        @test temp_vars.thermalcomfort isa ThermalComfort{FT,0}

        # Test field access for scalar case
        @test temp_vars.tempvec.TRoofImp === 0.0
        @test temp_vars.thermalcomfort.UTCI === 0.0
    end

    @testset "TimeSeries" begin
        temp_vars = initialize_temperature_variables(FT, TimeSeries(), hours)

        # Test structure
        @test temp_vars isa TemperatureVariables{FT,1}

        # Test components
        @test temp_vars.tempvec isa TempVec{FT,1}
        @test temp_vars.tempdamp isa TempDamp{FT,1}
        @test temp_vars.mrt isa MRT{FT,1}
        @test temp_vars.thermalcomfort isa ThermalComfort{FT,1}

        # Test field dimensions
        @test size(temp_vars.tempvec.TRoofImp) == (hours,)
        @test size(temp_vars.thermalcomfort.UTCI) == (hours,)
    end
end

@testset "get/setindex" begin
    hours = 24
    temp_vars = initialize_temperature_variables(FT, TimeSeries(), hours)

    TRoofImp = FT(300.0)
    TDampGroundImp = FT(290.0)
    LWR_Person = FT(310.0)
    UTCI = FT(25.0)

    x = temp_vars[1]
    x.tempvec.TRoofImp = TRoofImp
    x.tempdamp.TDampGroundImp = TDampGroundImp
    x.mrt.LWR_Person = LWR_Person
    x.thermalcomfort.UTCI = UTCI

    temp_vars[2] = x
    @test temp_vars.tempvec.TRoofImp[2] == TRoofImp
    @test temp_vars.tempdamp.TDampGroundImp[2] == TDampGroundImp
    @test temp_vars.mrt.LWR_Person[2] == LWR_Person
    @test temp_vars.thermalcomfort.UTCI[2] == UTCI
end
