using Test
using NCDatasets: NCDataset
using UrbanTethysChloris: create_model, initialize!, run_simulation, unpack_results
using UrbanTethysChloris: extended_outputs
using ..TestUtils: load_test_netcdf, load_test_parameters

FT = Float64
NN = 4

@testset "Output manager" begin
    @testset "no output file" begin
        model, forcing = create_model(
            FT,
            joinpath(@__DIR__, "data", "input_data.nc"),
            joinpath(@__DIR__, "data", "parameters.yaml"),
        )
        initialize!(model, forcing)
        manager, view_factor, view_factor_point = run_simulation(model, forcing; NN=NN)
        @test isnothing(manager)
    end

    @testset "streamed file" begin
        mktempdir() do tmp
            outpath = joinpath(tmp, "results.nc")
            model, forcing = create_model(
                FT,
                joinpath(@__DIR__, "data", "input_data.nc"),
                joinpath(@__DIR__, "data", "parameters.yaml"),
            )
            initialize!(model, forcing)
            manager, view_factor, view_factor_point = run_simulation(
                model,
                forcing;
                NN=NN,
                output_level=extended_outputs,
                output_filename=outpath,
            )
            @test !isnothing(manager)

            ds = NCDataset(outpath; maskingvalue=FT(NaN))
            @test haskey(ds.group, "temperature.tempvec")
            @test haskey(ds.group, "humidity.Results2m")
            @test haskey(ds.group, "waterflux.Owater")
            @test haskey(ds.group, "forcing.meteorological")
            @test haskey(ds.group, "parameters.urbangeometry")
            @test haskey(ds.group, "initial.OwaterInitial")
            @test ds.dim["rooflayers"] == model.parameters.soil.roof.ms
            @test ds.dim["groundlayers"] == model.parameters.soil.ground.ms

            # hourly values match the model state at the last saved step
            T2m = Float64.(coalesce.(Array(ds.group["humidity.Results2m"]["T2m"]), FT(NaN)))
            @test size(T2m, 1) >= NN
            @test T2m[NN] ≈ model.variables.humidity.Results2m.T2m

            Ow = Float64.(
                coalesce.(Array(ds.group["waterflux.Owater"]["OwGroundSoilVeg"]), FT(NaN))
            )
            @test size(Ow, 2) == model.parameters.soil.ground.ms
            @test Ow[NN, :] ≈ Array(model.variables.waterflux.Owater.OwGroundSoilVeg)

            # the static forcing series is written once at creation time
            Ta = Array(ds.group["forcing.meteorological"]["Tatm"])
            @test Ta ≈ Array(forcing.meteorological.Tatm)

            # the initial soil moisture is saved as a static variable
            OI = Array(ds.group["initial.OwaterInitial"]["OwGroundSoilVeg"])
            @test length(OI) == model.parameters.soil.ground.ms

            # the first-timestep energy balance is reset to zero, similar to MATLAB
            EB = Float64.(
                coalesce.(Array(ds.group["energybalance.EB"]["EBRoofImp"]), FT(NaN))
            )
            @test EB[1] == 0.0

            close(ds)

            # unpack_results reads the groups back with the former keys
            results = unpack_results(
                outpath, 1:NN, "temperature.tempvec", "humidity.Results2m"
            )
            @test results[:Results2m][:T2m] ≈ T2m[1:NN]
            @test results[:tempvec][:TRoofImp][NN] ≈
                model.variables.temperature.tempvec.TRoofImp
        end
    end
end
