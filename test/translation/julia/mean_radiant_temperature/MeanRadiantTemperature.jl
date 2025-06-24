using SafeTestsets

@safetestset "SWRDiffPerson" begin
    include("swr_diff_person.jl")
end

@safetestset "SWRDirPerson" begin
    include("swr_dir_person.jl")
end
