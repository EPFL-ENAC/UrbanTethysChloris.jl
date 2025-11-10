using SafeTestsets

@safetestset "ForcingInputs" begin
    include(joinpath("forcinginputs", "ForcingInputs.jl"))
end

@safetestset "Parameters" begin
    include(joinpath("parameters", "Parameters.jl"))
end

@safetestset "ModelVariables" begin
    include(joinpath("modelvariables", "ModelVariables.jl"))
end
