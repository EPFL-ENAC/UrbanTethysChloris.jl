using SafeTestsets

@safetestset "SWRDiffPerson" begin
    include("swr_diff_person.jl")
end

@safetestset "SWRDirPerson" begin
    include("swr_dir_person.jl")
end

@safetestset "PersonInShade" begin
    include("person_in_shade.jl")
end

@safetestset "MeanRadiantTemperature" begin
    include("mean_radiant_temperature.jl")
end
