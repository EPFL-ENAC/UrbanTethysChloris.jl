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
    "Radius_tree" => 0.289,
    "Hcan_max" => NaN,
    "Hcan_std" => NaN,
    "trees" => true,
    "ftree" => 1.0,
)

data["urbangeometry"]["Height_tree"] = 5.0 - data["urbangeometry"]["Radius_tree"]
data["urbangeometry"]["Distance_tree"] = 0.5 + data["urbangeometry"]["Radius_tree"]

data["surfacefractions"] = Dict{String,Any}()
data["surfacefractions"]["roof"] = Dict{String,Any}(
    "fveg" => 0.5, "fimp" => 0.5, "Per_runoff" => 1.0
)
data["surfacefractions"]["ground"] = Dict{String,Any}(
    "fveg" => 0.545, "fbare" => 0.0, "fimp" => 0.455, "Per_runoff" => 0.9
)

# Vegetation
data["vegetation"] = Dict{String,Any}()
## Roof
data["vegetation"]["roof"] = Dict{String,Any}(
    "LAI" => 2.5,
    "SAI" => 0.001,
    "hc" => 0.15,
    "d_leaf" => 0.8,
    "CASE_ROOT" => 1,
    "ZR95" => 70.0,
    "ZR50" => NaN,
    "ZRmax" => NaN,
    "Rrootl" => 3800.0,
    "PsiL50" => -4.0,
    "PsiX50" => -4.5,
    "FI" => 0.081,
    "Do" => 1000.0,
    "a1" => 6.0,
    "go" => 0.01,
    "CT" => 3,
    "DSE" => 0.656,
    "Ha" => 55.0,
    "gmes" => Inf,
    "rjv" => 2.4,
    "Kopt" => 0.5,
    "Knit" => 0.15,
    "Vmax" => 68.0,
    "mSl" => 0.0,
    "e_rel" => 1.0,
    "e_relN" => 1.0,
    "Psi_sto_00" => -0.5,
    "Psi_sto_50" => -3.0,
    "Sl" => 0.035,
)

data["vegetation"]["roof"]["h_disp"] = 2.0 / 3.0 * data["vegetation"]["roof"]["hc"]

## Ground
data["vegetation"]["ground"] = data["vegetation"]["roof"]
data["vegetation"]["ground"]["ZR95"] = 250.0

## Tree
data["vegetation"]["tree"] = data["vegetation"]["roof"]
data["vegetation"]["tree"]["LAI"] = 5.0
data["vegetation"]["tree"]["SAI"] = 0.2
data["vegetation"]["tree"]["hc"] = NaN
data["vegetation"]["tree"]["ZR95"] = 1000.0
data["vegetation"]["tree"]["Rrootl"] = 4000.0
data["vegetation"]["tree"]["PsiL50"] = -3.0
data["vegetation"]["tree"]["a1"] = 9.0
data["vegetation"]["tree"]["DSE"] = 0.649
data["vegetation"]["tree"]["Ha"] = 76.0
data["vegetation"]["tree"]["Sl"] = 0.024
data["vegetation"]["tree"]["SPARTREE"] = 2

# Thermal properties
data["thermal"] = Dict{String,Any}()
data["thermal"]["roof"] = Dict{String,Any}("lan_dry" => 0.67, "cv_s" => 1e6)

data["thermal"]["ground"] = Dict{String,Any}("lan_dry" => 1.2, "cv_s" => 1.5e6)

data["thermal"]["wall"] = data["thermal"]["roof"]

YAML.write_file(joinpath(@__DIR__, "..", "data", "parameters.yaml"), data)
YAML.write_file(joinpath(@__DIR__, "..", "test", "data", "parameters.yaml"), data)
