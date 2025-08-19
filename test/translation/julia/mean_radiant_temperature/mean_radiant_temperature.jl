using Test
using MAT
using UrbanTethysChloris.MeanRadiantTemperature: mean_radiant_temperature
using UrbanTethysChloris.RayTracing: ViewFactorPoint
using ....TestUtils:
    create_urban_geometry_parameters,
    create_height_dependent_vegetation_parameters,
    load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("MRT.MeanRadiantTemperature.mat")

# Create parameter structs
ParVegTree = create_height_dependent_vegetation_parameters(
    FT; LAI=input_vars["ParVegTree"]["LAI"], Kopt=input_vars["ParVegTree"]["Kopt"]
)

Geometry_m = create_urban_geometry_parameters(
    FT;
    trees=Bool(input_vars["ParTree"]["trees"]),
    hcanyon=input_vars["geometry"]["hcanyon"],
    wcanyon=input_vars["geometry"]["wcanyon"],
    htree=input_vars["geometry"]["htree"],
    radius_tree=input_vars["geometry"]["radius_tree"],
    distance_tree=input_vars["geometry"]["distance_tree"],
    Width_canyon=input_vars["Gemeotry_m"]["Width_canyon"],
)

SunPosition = (;
    theta_Z=input_vars["SunPosition"]["theta_Z"],
    theta_n=input_vars["SunPosition"]["theta_n"],
    zeta_S=input_vars["SunPosition"]["zeta_S"],
    TimeOfMaxSolAlt=input_vars["SunPosition"]["TimeOfMaxSolAlt"],
    Datam=vec(input_vars["SunPosition"]["Datam"]),
)

Person = (;
    PositionPz=input_vars["Person"]["PositionPz"],
    PositionPx=input_vars["Person"]["PositionPx"],
)

SWRout_t = (;
    SWRoutTotalGround=input_vars["SWRout_t"]["SWRoutTotalGround"],
    SWRoutTree=input_vars["SWRout_t"]["SWRoutTree"],
    SWRoutWallSun=input_vars["SWRout_t"]["SWRoutWallSun"],
    SWRoutWallShade=input_vars["SWRout_t"]["SWRoutWallShade"],
)

LWRout_t = (;
    LWRoutTotalGround=input_vars["LWRout_t"]["LWRoutTotalGround"],
    LWRoutTree=input_vars["LWRout_t"]["LWRoutTree"],
    LWRoutWallSun=input_vars["LWRout_t"]["LWRoutWallSun"],
    LWRoutWallShade=input_vars["LWRout_t"]["LWRoutWallShade"],
)

MeteoData = (;
    SW_dir=input_vars["MeteoData"]["SW_dir"],
    SW_diff=input_vars["MeteoData"]["SW_diff"],
    LWR=input_vars["MeteoData"]["LWR"],
    SunDSM_MRT=input_vars["MeteoData"]["SunDSM_MRT"],
)

vfp = ViewFactorPoint{FT}(;
    F_pg=input_vars["ViewFactorPoint"]["F_pg"],
    F_ps=input_vars["ViewFactorPoint"]["F_ps"],
    F_pt=input_vars["ViewFactorPoint"]["F_pt"],
    F_pwLeft=input_vars["ViewFactorPoint"]["F_pwLeft"],
    F_pwRight=input_vars["ViewFactorPoint"]["F_pwRight"],
)

@testset "MATLAB" begin
    Tmrt, BooleanInSun, SWRdir_Person, SWRdir_in_top, SWRdir_in_bottom, SWRdir_in_east, SWRdir_in_south, SWRdir_in_west, SWRdir_in_north, SWRdiff_Person, LWR_Person = mean_radiant_temperature(
        SWRout_t, LWRout_t, MeteoData, vfp, ParVegTree, Geometry_m, SunPosition, Person
    )

    @test Tmrt ≈ output_vars["Tmrt"]
    @test BooleanInSun ≈ output_vars["BoleanInSun"]
    @test SWRdir_Person ≈ output_vars["SWRdir_Person"]
    @test SWRdir_in_top ≈ output_vars["SWRdir_in_top"]
    @test SWRdir_in_bottom ≈ output_vars["SWRdir_in_bottom"]
    @test SWRdir_in_east ≈ output_vars["SWRdir_in_east"]
    @test SWRdir_in_south ≈ output_vars["SWRdir_in_south"]
    @test SWRdir_in_west ≈ output_vars["SWRdir_in_west"]
    @test SWRdir_in_north ≈ output_vars["SWRdir_in_north"]
    @test SWRdiff_Person ≈ output_vars["SWRdiff_Person"]
    @test LWR_Person ≈ output_vars["LWR_Person"]
end
