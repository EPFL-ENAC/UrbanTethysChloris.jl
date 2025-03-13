using NCDatasets
using YAML
using MAT

FT = Float64

# Check if data directory exists, if not create it
data_dir = joinpath(@__DIR__, "..", "data")
!isdir(data_dir) && mkdir(data_dir)

# Define files and their GitHub URLs
repo_url = "https://github.com/NaikaMeili/UTC_ModelCode/raw/b7c4c0617133681b678ec066cbcd881a8fb97aae/UTC_Model/"
files = Dict(
    "ForcingData_ZH2010.mat" => repo_url * "%2Bdata_functions/ForcingData_ZH2010.mat"
)

# Check each file and download if missing
for (file, url) in files
    filepath = joinpath(data_dir, file)
    if !isfile(filepath)
        @info "Downloading $file..."
        download(url, filepath)
    end
end

# Create the parameters file
data = Dict{String,Any}()

data["urbangeometry"] = Dict{String,Any}(
    "Height_canyon" => 6.5,
    "Width_canyon" => 5.78,
    "Width_roof" => 4.73,
    "Radius_tree" => 100.0,
    "Height_tree" => 0.289,
    "Distance_tree" => 100.0,
    "Hcan_max" => NaN,
    "Hcan_std" => NaN,
    "trees" => true,
    "ftree" => 1.0,
)

data["urbangeometry"]["Height_tree"] = 5.0 - data["urbangeometry"]["Radius_tree"]
data["urbangeometry"]["Distance_tree"] = 0.5 + data["urbangeometry"]["Radius_tree"]

data["surfacefractions"] = Dict{String,Any}(
    "fveg_R" => 0.5,
    "fimp_R" => 0.5,
    "Per_runoff_R" => 1.0,
    "fveg_G" => 0.545,
    "fbare_G" => 0.0,
    "fimp_G" => 0.455,
    "Per_runoff_G" => 0.9,
)

YAML.write_file(joinpath(@__DIR__, "..", "data", "parameters.yaml"), data)
