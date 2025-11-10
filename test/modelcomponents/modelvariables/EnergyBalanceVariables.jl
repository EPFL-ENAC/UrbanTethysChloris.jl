using Test
using UrbanTethysChloris.ModelComponents.ModelVariables:
    WBRoof,
    initialize_wbroof,
    WBCanyonIndv,
    initialize_wbcanyon_indv,
    WBCanyonTot,
    initialize_wbcanyon_tot,
    EB,
    initialize_eb,
    SolverVariables,
    initialize_solver_variables,
    EnergyBalanceVariables,
    initialize_energy_balance_variables,
    TimeSlice,
    TimeSeries

FT = Float64

@testset "Subsets" begin
    @testset "WBRoof" begin
        @testset "TimeSlice" begin
            wbroof = initialize_wbroof(FT, TimeSlice())
            @test wbroof isa WBRoof{FT,0}
            @test wbroof.WBRoofImp === 0.0
            @test wbroof.WBRoofVeg === 0.0
            @test wbroof.WBRoofTot === 0.0
            for field in fieldnames(WBRoof)
                @test isa(getproperty(wbroof, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            wbroof = initialize_wbroof(FT, TimeSeries(), hours)
            @test wbroof isa WBRoof{FT,1}
            for field in fieldnames(WBRoof)
                @test isa(getproperty(wbroof, field), Array{FT,1})
                @test size(getproperty(wbroof, field)) == (hours,)
                @test getproperty(wbroof, field)[1] == 0
            end
        end
    end

    @testset "WBCanyonIndv" begin
        @testset "TimeSlice" begin
            wbcanyon_indv = initialize_wbcanyon_indv(FT, TimeSlice())
            @test wbcanyon_indv isa WBCanyonIndv{FT,0}
            @test wbcanyon_indv.WB_In_tree === 0.0
            @test wbcanyon_indv.WB_In_gveg === 0.0
            @test wbcanyon_indv.WB_Soil_gveg === 0.0
            for field in fieldnames(WBCanyonIndv)
                @test isa(getproperty(wbcanyon_indv, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            wbcanyon_indv = initialize_wbcanyon_indv(FT, TimeSeries(), hours)
            @test wbcanyon_indv isa WBCanyonIndv{FT,1}
            for field in fieldnames(WBCanyonIndv)
                @test isa(getproperty(wbcanyon_indv, field), Array{FT,1})
                @test size(getproperty(wbcanyon_indv, field)) == (hours,)
                @test getproperty(wbcanyon_indv, field)[1] == 0
            end
        end
    end

    @testset "WBCanyonTot" begin
        @testset "TimeSlice" begin
            wbcanyon_tot = initialize_wbcanyon_tot(FT, TimeSlice())
            @test wbcanyon_tot isa WBCanyonTot{FT,0}
            @test wbcanyon_tot.WBsurf_tree === 0.0
            @test wbcanyon_tot.WBsoil_veg === 0.0
            @test wbcanyon_tot.WBcanyon_level === 0.0
            for field in fieldnames(WBCanyonTot)
                @test isa(getproperty(wbcanyon_tot, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            wbcanyon_tot = initialize_wbcanyon_tot(FT, TimeSeries(), hours)
            @test wbcanyon_tot isa WBCanyonTot{FT,1}
            for field in fieldnames(WBCanyonTot)
                @test isa(getproperty(wbcanyon_tot, field), Array{FT,1})
                @test size(getproperty(wbcanyon_tot, field)) == (hours,)
                @test getproperty(wbcanyon_tot, field)[1] == 0
            end
        end
    end

    @testset "EB" begin
        @testset "TimeSlice" begin
            eb = initialize_eb(FT, TimeSlice())
            @test eb isa EB{FT,0}
            @test eb.EBRoofImp === 0.0
            @test eb.EBGroundVeg === 0.0
            @test eb.EBCanyonQ === 0.0
            for field in fieldnames(EB)
                @test isa(getproperty(eb, field), FT)
            end
        end
        @testset "TimeSeries" begin
            hours = 24
            eb = initialize_eb(FT, TimeSeries(), hours)
            @test eb isa EB{FT,1}
            for field in fieldnames(EB)
                @test isa(getproperty(eb, field), Array{FT,1})
                @test size(getproperty(eb, field)) == (hours,)
                @test getproperty(eb, field)[1] == 0
            end
        end
    end

    @testset "SolverVariables" begin
        @testset "TimeSlice" begin
            sv = initialize_solver_variables(FT, TimeSlice())
            @test sv isa SolverVariables{FT,0,1}
            @test sv.Success === false
            @test all(sv.ValuesEB .== 0.0)
            @test all(sv.Tsolver .== 0.0)
            for field in setdiff(fieldnames(SolverVariables), [:Success])
                @test isa(getproperty(sv, field), Vector{FT})
            end
            @test isa(getproperty(sv, :Success), Bool)
        end
        @testset "TimeSeries" begin
            hours = 24
            sv = initialize_solver_variables(FT, TimeSeries(), hours)
            @test sv isa SolverVariables{FT,1,2}
            for field in setdiff(fieldnames(SolverVariables), [:Success])
                @test isa(getproperty(sv, field), Matrix{FT})
                @test size(getproperty(sv, field)) == (hours, 22)
                @test all(getproperty(sv, field)[1, :] .== 0)
            end
            @test isa(getproperty(sv, :Success), Vector{Bool})
            @test size(getproperty(sv, :Success)) == (hours,)
            @test all(getproperty(sv, :Success) .== false)
        end
    end
end

@testset "EnergyBalanceVariables" begin
    @testset "TimeSlice" begin
        energy_balance_vars = initialize_energy_balance_variables(FT, TimeSlice())

        # Test structure
        @test energy_balance_vars isa EnergyBalanceVariables{FT,0}

        # Test fields are properly initialized
        @test isa(energy_balance_vars.WBRoof, WBRoof{FT,0})
        @test isa(energy_balance_vars.WBCanyonIndv, WBCanyonIndv{FT,0})
        @test isa(energy_balance_vars.WBCanyonTot, WBCanyonTot{FT,0})
        @test isa(energy_balance_vars.EB, EB{FT,0})
        @test isa(energy_balance_vars.Solver, SolverVariables{FT,0,1})

        # Test some random fields
        @test energy_balance_vars.WBRoof.WBRoofImp === 0.0
        @test energy_balance_vars.WBCanyonIndv.WB_In_tree === 0.0
        @test energy_balance_vars.WBCanyonTot.WBsoil_veg === 0.0
        @test energy_balance_vars.EB.EBRoofImp === 0.0
        @test energy_balance_vars.Solver.Success === false
    end

    @testset "TimeSeries" begin
        hours = 24
        energy_balance_vars = initialize_energy_balance_variables(FT, TimeSeries(), hours)

        # Test structure
        @test energy_balance_vars isa EnergyBalanceVariables{FT,1,2}

        # Test fields are properly initialized
        @test isa(energy_balance_vars.WBRoof, WBRoof{FT,1})
        @test isa(energy_balance_vars.WBCanyonIndv, WBCanyonIndv{FT,1})
        @test isa(energy_balance_vars.WBCanyonTot, WBCanyonTot{FT,1})
        @test isa(energy_balance_vars.EB, EB{FT,1})
        @test isa(energy_balance_vars.Solver, SolverVariables{FT,1,2})

        # Test dimensions of a few fields
        @test size(energy_balance_vars.WBRoof.WBRoofImp) == (hours,)
        @test size(energy_balance_vars.WBCanyonIndv.WB_In_tree) == (hours,)
        @test size(energy_balance_vars.WBCanyonTot.WBsoil_veg) == (hours,)
        @test size(energy_balance_vars.EB.EBRoofImp) == (hours,)
        @test size(energy_balance_vars.Solver.Success) == (hours,)

        # Test initialization values
        @test energy_balance_vars.WBRoof.WBRoofImp[1] == 0.0
        @test energy_balance_vars.WBCanyonIndv.WB_In_tree[1] == 0.0
        @test energy_balance_vars.WBCanyonTot.WBsoil_veg[1] == 0.0
        @test energy_balance_vars.EB.EBRoofImp[1] == 0.0
        @test energy_balance_vars.Solver.Success[1] == false
    end
end

@testset "get/setindex" begin
    hours = 24
    energy_balance_vars = initialize_energy_balance_variables(FT, TimeSeries(), hours)

    WBRoofImp = FT(1.5)
    WB_In_tree = FT(3.5)
    WBsoil_veg = FT(4.5)
    EBRoofImp = FT(2.5)
    ValuesEB = fill(FT(1.5), 22)

    # Get index 1
    ebv = energy_balance_vars[1]
    ebv.WBRoof.WBRoofImp = WBRoofImp
    ebv.WBCanyonIndv.WB_In_tree = WB_In_tree
    ebv.WBCanyonTot.WBsoil_veg = WBsoil_veg
    ebv.EB.EBRoofImp = EBRoofImp
    ebv.Solver.ValuesEB = ValuesEB
    energy_balance_vars[2] = ebv

    @test energy_balance_vars.WBRoof.WBRoofImp[2] == WBRoofImp
    @test energy_balance_vars.WBCanyonIndv.WB_In_tree[2] == WB_In_tree
    @test energy_balance_vars.WBCanyonTot.WBsoil_veg[2] == WBsoil_veg
    @test energy_balance_vars.EB.EBRoofImp[2] == EBRoofImp
    @test energy_balance_vars.Solver.ValuesEB[2, :] == ValuesEB
end
