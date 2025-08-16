using Test
using MAT
using UrbanTethysChloris.MeanRadiantTemperature: swr_diff_person
using UrbanTethysChloris.RayTracing: ViewFactorPoint
using ....TestUtils: load_matlab_data

FT = Float64
input_vars, output_vars = load_matlab_data("MRT.SWRDiffPerson.mat")

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
    SW_diff=input_vars["MeteoData"]["SW_diff"], LWR=input_vars["MeteoData"]["LWR"]
)
vfp = ViewFactorPoint{FT}(;
    F_pg=input_vars["ViewFactorPoint"]["F_pg"],
    F_ps=input_vars["ViewFactorPoint"]["F_ps"],
    F_pt=input_vars["ViewFactorPoint"]["F_pt"],
    F_pwLeft=input_vars["ViewFactorPoint"]["F_pwLeft"],
    F_pwRight=input_vars["ViewFactorPoint"]["F_pwRight"],
)
TimeOfMaxSolAlt = input_vars["TimeOfMaxSolAlt"]
TimeHr = input_vars["TimeHr"]

@testset "MATLAB" begin
    SWRdiff_Person, LWR_Person = swr_diff_person(
        SWRout_t, LWRout_t, MeteoData, vfp, TimeOfMaxSolAlt, TimeHr
    )
    @test SWRdiff_Person == output_vars["SWRdiff_Person"]
    @test LWR_Person == output_vars["LWR_Person"]
end
