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
data["thermal"]["tree"] = Dict{String,Any}("Cthermal_leaf" => 640.0)

data["optical"] = Dict{String,Any}()
data["optical"]["wall"] = Dict{String,Any}("albedo" => 0.4, "emissivity" => 0.95)
data["optical"]["tree"] = Dict{String,Any}("albedo" => 0.2, "emissivity" => 0.95)

LAI_R = data["vegetation"]["roof"]["LAI"]
SAI_R = data["vegetation"]["roof"]["SAI"]
data["optical"]["roof"] = Dict{String,Any}(
    "aveg" => 0.2, "aimp" => 0.15, "eveg" => 1 - exp(-(LAI_R + SAI_R)), "eimp" => 0.95
)

LAI_G = data["vegetation"]["ground"]["LAI"]
SAI_G = data["vegetation"]["ground"]["SAI"]
data["optical"]["ground"] = Dict{String,Any}(
    "aveg" => 0.2,
    "abare" => 0.15,
    "aimp" => 0.1,
    "eveg" => 1 - exp(-(LAI_G + SAI_G)),
    "ebare" => 0.95,
    "eimp" => 0.95,
)

data["soil"] = Dict{String,Any}()

# Roof soil parameters
data["soil"]["roof"] = Dict{String,Any}(
    "Pcla" => 0.20,
    "Psan" => 0.40,
    "Porg" => 0.025,
    "In_max_imp" => 0.25,
    "In_max_ground" => 10.0,
    "Sp_In" => 0.2,
    "Kimp" => 0.0,
    "Kfc" => 0.2,
    "Phy" => 10000.0,
    "SPAR" => 2,
    "Kbot" => NaN,
    "dz1" => 0.1,
    "dz2" => 0.1,
    "Zs" => [0.0, 10.0, 20.0, 50.0, 100.0],
    "FixSM" => true,
    "FixSM_LayerStart" => 1,
    "FixSM_LayerEnd" => 4,
)

# Ground soil parameters
data["soil"]["ground"] = Dict{String,Any}(
    "Pcla" => 0.20,
    "Psan" => 0.40,
    "Porg" => 0.025,
    "In_max_imp" => 0.5,
    "In_max_underveg" => 10.0,
    "In_max_bare" => 10.0,
    "Sp_In" => 0.2,
    "Kimp" => 0.001,
    "Kfc" => 0.2,
    "Phy" => 10000.0,
    "SPAR" => 2,
    "Kbot" => NaN,
    "Zs" => [
        0.0,
        10.0,
        20.0,
        50.0,
        100.0,
        150.0,
        200.0,
        300.0,
        400.0,
        600.0,
        800.0,
        1000.0,
        1500.0,
        2000.0,
    ],
    "FixSM" => true,
    "FixSM_LayerStart" => 6,
    "FixSM_LayerEnd" => 13,
)

# Wall soil parameters
data["soil"]["wall"] = Dict{String,Any}("dz1" => 0.1, "dz2" => 0.1)

# Tree interception parameter
data["soil"]["Sp_In_T"] = 0.2

data["building_energy"] = Dict{String,Any}()
data["building_energy"]["indoor_optical"] = Dict{String,Any}(
    "abc" => 0.3,
    "abw" => 0.3,
    "abg" => 0.3,
    "abm" => 0.3,
    "ec" => 0.95,
    "eg" => 0.95,
    "ew" => 0.95,
    "em" => 0.95,
)
data["building_energy"]["thermal"] = Dict{String,Any}(
    "IntMassOn" => 0,
    "FloorHeight" => 3.0,
    "dzFloor" => 0.2,
    "dzWall" => 0.2,
    "lan_ground_floor" => 1.2,
    "cv_ground_floor" => 1.5e6,
    "lan_floor_IntMass" => 0.67,
    "cv_floor_IntMass" => 1.0e6,
    "lan_wall_IntMass" => 0.67,
    "cv_wall_IntMass" => 1.0e6,
)
data["building_energy"]["windows"] = Dict{String,Any}(
    "WindowsOn" => 1,
    "GlazingRatio" => 0.15,
    "Uvalue" => 4.95,
    "lan_windows" => NaN,
    "cv_glass" => 2.1e6,
    "dztot" => 0.02,
    "SHGC" => 0.8,
    "SolarTransmittance" => 0.6,
    "SolarAbsorptivity" => 0.0,
    "SolarAlbedo" => 0.4,
)
data["building_energy"]["hvac"] = Dict{String,Any}(
    "ACon" => 1,
    "Heatingon" => 1,
    "TsetpointCooling" => 298.15,
    "TsetpointHeating" => 293.15,
    "RHsetpointCooling" => 60.0,
    "RHsetpointHeating" => NaN,
    "ACH" => 0.5,
    "COPAC" => 3.26,
    "COPHeat" => 0.9,
    "f_ACLatentToQ" => 1.0,
)

# Add person parameters
data["person"] = Dict{String,Any}(
    "PositionPx" => data["urbangeometry"]["Width_canyon"] / 2,
    "PositionPz" => 1.1,
    "PersonWidth" => 0.06 / 2,
    "PersonHeight" => 0.22 / 2,
    "HeightWind" => 1.1,
)

YAML.write_file(joinpath(@__DIR__, "..", "data", "parameters.yaml"), data)
YAML.write_file(joinpath(@__DIR__, "..", "test", "data", "parameters.yaml"), data)
